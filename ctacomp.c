#include <stdio.h>
#include <stdlib.h>
#include "fitsio2.h"
#include <string.h>

//Structure helper for blocks headers and compresion schemes
struct BlockHeader
{
    unsigned long long size;
    char          ordering;
    unsigned char numProcs;
    unsigned short      processings[];
} __attribute__((__packed__));

//struct LUT;

struct LUT
{
    unsigned short symbol;
    unsigned char  nbits;
    int            is_leaf;
    struct LUT*    lut;
} LUT;

void erase_lut(struct LUT* lut)
{
    int i;
    for (i=0;i<256;i++)
        if (lut[i].lut != NULL)
            erase_lut(lut[i].lut);
    free(lut);
}


void copy_swap(char *dest, const char *src, int num, int elem_size)
{
    int i;
    const char *pend = src + num*elem_size;
    for (const char *ptr = src; ptr<pend; ptr+=elem_size, dest+=elem_size)
        for (i=0;i<elem_size;i++)
            dest[i] = ptr[elem_size-i];
}


void add_symbol(struct LUT** lut, unsigned short sym, unsigned char nbits, size_t bits)
{
    int i;
    // If not lut available, create one
    if (!*lut) {
        (*lut) = malloc(256 * sizeof(struct LUT));
        for (i=0;i<256;i++) {
            (*lut)[i].is_leaf = 0;
            (*lut)[i].lut     = NULL;
            (*lut)[i].symbol  = 0;
        } 
    }
    if (nbits>8) {
        add_symbol(&((*lut)[bits&0xff].lut), sym, nbits-8, bits>>8);
        return;
    }
    
    // From here, we deal with all the 8-bits leaves that contain the remaining code
    int n_leaves = 1<<(8-nbits);
    for (i=0;i<n_leaves;i++) {
        const unsigned char key = bits | (i<<nbits);
        (*lut)[key].symbol = sym;
        (*lut)[key].is_leaf = 1;
        (*lut)[key].nbits   = nbits;
    }
}

#define FACT_HUFFMAN_16 0x2
#define CTA_DIFF_16     0x3

int fits_ctadecomp(unsigned char* c, // Input buffer
                   unsigned long  input_len, //length of input in bytes
                   unsigned char array[], // output array
                   unsigned long  output_size, // size of the output array in bytes
                   int col_type, 
                   int col_width)
{
    
    printf("Dealing with column type %d of width %d\n", col_type, col_width);
    unsigned long output_num = output_size / 2;
    unsigned char* max_input_address = c + input_len;

    struct BlockHeader* head = (struct BlockHeader*)(c);

    int i;

    c += sizeof(struct BlockHeader) + head->numProcs*sizeof(unsigned short);

    if (head->numProcs != 2)
    {
        printf("ERROR: Number of processings applied to CTA compression not 2 as expected.");
        ffpmsg("ERROR: Number of processings applied to CTA compression not 2 as expected.");
        return 1;
    }
    
    if (head->processings[0] != CTA_DIFF_16 || head->processings[1] != FACT_HUFFMAN_16)
    {
        ffpmsg("ERROR: CTA compression processings not 0x3 and 0x2 as expected.");
        printf("ERROR: CTA compression processings not 0x3 and 0x2 as expected.");
        return 2;
    }

    // Read the size of the compressed data
    unsigned int buffer_length = 0;
    memcpy(&buffer_length, c, sizeof(buffer_length));
    c += sizeof(buffer_length);

    // Read the number of 16 bits data this encoding represents.
    unsigned long long data_count = 0;
    memcpy(&data_count, c, sizeof(data_count));
    c += sizeof(data_count);

    // Read the number of symbols in the look-up table
    unsigned long long sym_count = 0;
    memcpy(&sym_count, c, sizeof(sym_count));
    c += sizeof(sym_count);


    int total_bytes_read = 0;
    unsigned char* startup_c = c;
    // Read the look up table of the encoded symbols
    struct LUT* lut = NULL;
    for (i=0;i<sym_count;i++) {
        unsigned short symbol;
        memcpy(&symbol, c, sizeof(unsigned short));
        c += sizeof(unsigned short);
        if (sym_count == 1) {
            add_symbol(&lut, symbol, 0, 0);
            break;
        }
        unsigned char num_bits;
        memcpy(&num_bits, c, sizeof(unsigned char));
        c += sizeof(unsigned char);
        const unsigned char num_bytes = num_bits / 8 + (num_bits % 8 ? 1 : 0);
        if (num_bytes > sizeof(size_t)) {
            printf(">>>PROBLEM HERE %d %d %d i=%d\n\n", num_bytes, num_bits, sizeof(size_t), i);
            ffpmsg("ERROR: Huffman encoding length too large.");
            exit(-1);
            return 1;
        }
        //printf("Reading symbol #%d of %d bits on %d bytes (%d bytes read)\n", i, num_bits, num_bytes, total_bytes_read);
        total_bytes_read += num_bytes;
        size_t bits = 0;
        memcpy(&bits, c, num_bytes);
        c += num_bytes;

        add_symbol(&lut, symbol, num_bits, bits);
    }

    //if we have only one input symbol, fill up the output with it
    if (sym_count == 1) {
        for (i=0;i<output_num;i++) {
            ((unsigned short*)(array))[i] = lut->symbol;
        }
        return 0;
    }
   
    //take a temporary array to reshuffle the uncompressed data
    unsigned short* tmp_array = malloc(data_count*sizeof(unsigned short));


    //otherwise, until we exhaust the input bits
    struct LUT top_lut;
    top_lut.lut = lut;
    top_lut.is_leaf = 0;
    struct LUT* clut = &top_lut; //lut;
    unsigned char curbit = 0;
    unsigned long num_written = 0;
    while (c < max_input_address  && 
           num_written < data_count) {
        
        //look at the current input byte. Take a two-bytes word as input to make sure that we can extract a full byte
        const unsigned short *two = (unsigned short*)c;
        const unsigned char curbyte =  (*two >> curbit);
        
        //if we end up nowhere, there is a problem with the encoding
        if (!clut) {
            ffpmsg("Ended up with a NULL lut.");
            printf("ERROR: Ended up with a NULL lut. \n");
            return 1;
        }

        //otherwise, traverse the graph of symbols following the new input set of bits (curbyte)
        clut = &(clut->lut[curbyte]);

        //if we are not at a leaf yet, continue to traverse
        if (!clut->is_leaf)
        {
            c++;
            continue;
        }

        //if we've hit a leaf, assign the symbol to the output, and consume the corresponding bits
        tmp_array[num_written++] = clut->symbol;

        curbit    += clut->nbits;
        //and start over from the root of the graph (lut)
        clut = &top_lut;
        if (curbit>=8){
            curbit %= 8;
            c++;
        }
    }

    printf("Writing done. num written: %lu, output_num: %lu, data_count: %u\n", num_written, output_num, data_count);
    printf("\n\n");

    // Un-apply the subtraction of the previous element
    for (i=1;i<data_count;i++)
        tmp_array[i] += tmp_array[i-1];

    // Copy the temp array data to the output array
    if (col_width == 1) {
        if (data_count != output_num) {
            printf("PROBLEM HERE!!! NON-ARRAY SIZE DISAGREE %d vs %d\n", data_count , output_num);
            exit(-1);
        }
        memcpy(array, tmp_array, output_size);
    }
    else { /* Deal with arrays of possibly variable length*/
        int src_index = 0;
        int dst_index = 0;
        unsigned char* char_array = (unsigned char*)(tmp_array);
        while (src_index < (data_count*sizeof(unsigned short)-4)) { /* -4 to make sure to at least be able to read the expected next size*/ 
            // FIXME What the #$%$## is this -4 really ? Numbers should match !!!
            int this_row_bytes = 0;
            memcpy(&this_row_bytes, &char_array[src_index], sizeof(this_row_bytes));
            src_index += sizeof(this_row_bytes);
            int bytes_factor = 1;
            switch (col_type) {
                case TBIT:
                    //FIXME Figure out what to do with TBITS
                    printf("Not sure what to do with TBITS...\n");
                    exit(-1);
                break;
                case TBYTE:
                case TSBYTE:
                case TSTRING:
                break;
                case TUSHORT:
                case TSHORT:
                    bytes_factor = 2;
                break;
                case TUINT:
                case TINT:
                case TFLOAT:
                    bytes_factor = 4;
                break;
                case TULONG:
                case TLONG:
                case TULONGLONG:
                case TLONGLONG:
                case TDOUBLE:
                case TCOMPLEX:
                case TDBLCOMPLEX:
                    bytes_factor = 8;
                break;
                default:
                    printf("ERROR: Got a default case in ctacomp.c. col_type was %d\n", col_type);
                    exit(-1);
                break;
            };
            printf("Copying from %d to %d size %d\n", src_index, dst_index, this_row_bytes);
            copy_swap(&(array[dst_index]), &(char_array[src_index]), this_row_bytes / bytes_factor, bytes_factor);
            dst_index += this_row_bytes;
            src_index += this_row_bytes;
            // Pad any row shorter than the expected size with zeros
            for (i=this_row_bytes;i<col_width*bytes_factor;i++)
                array[dst_index++] = 0;
        }
        printf("Reduced %d elems from tmp_array to bytes array of size %d\n", src_index, dst_index);
        printf("Expected sizes were %d and %d\n", data_count*2, output_size);
    }



    free(tmp_array);

    // swap the bytes if the ordering was big endian


    erase_lut(lut);

    return 0;
}


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
                   unsigned short array[], // output array
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
            array[i] = lut->symbol;
        }
        return 0;
    }
    int j;
    int bytes_to_go = max_input_address - c;

    //otherwise, until we exhaust the input bits
    struct LUT top_lut;
    top_lut.lut = lut;
    top_lut.is_leaf = 0;
    struct LUT* clut = &top_lut; //lut;
    unsigned char curbit = 0;
    unsigned long num_written = 0;
    while (num_written < output_num && 
           c < max_input_address  && 
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
        array[num_written++] = clut->symbol;

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
    for (i=1;i<output_num;i++)
        array[i] += array[i-1];

    erase_lut(lut);

    return 0;
}


/*
 Original Algorithm by Thomas Bretz and Etienne Lyard for the FACT telescope.
 Now available at the CTA gitlab: 
 https://gitlab.cta-observatory.org/cta-computing/common/acada-array-elements/adh-apis/-/blob/master/zfits/FlatProtobufZOFits.cpp?ref_type=heads#L2411
 
 Recoding in C by Etienne Lyard etienne.lyard@unige.ch
*/

#include <stdio.h>
#include <stdlib.h>
#include "fitsio2.h"
#include <string.h>

/* Structure helper for compressed blocks headers */
struct BlockHeader
{
    unsigned long long size;
    char          ordering;
    unsigned char numProcs;
    unsigned short      processings[];
} __attribute__((__packed__));

/* Structure holding the huffman encoding */
struct LUT
{
    unsigned short symbol;
    unsigned char  nbits;
    int            is_leaf;
    struct LUT*    lut;
} LUT;

/* Erase the huffman coding recursively
   @param lut the top-level look-up-table node */
void erase_lut(struct LUT* lut)
{
    int i;
    for (i=0;i<256;i++)
        if (lut[i].lut != NULL)
            erase_lut(lut[i].lut);
    free(lut);
}

/* Apply any offset defined in the FITS header
   @param array     The memory containing the values to offset
   @param num       The number of values contained in the array
   @param elem_size The size in bytes of a given value
   @param offset    The value to use as offset */
void apply_offset(char* array, 
                  int   num, 
                  int   elem_size, 
                  unsigned long long offset)
{
    if (offset == 0)
        return;

    int i;
    switch (elem_size)
    {
        case 1: ;
            unsigned char* uchar_array = (unsigned char*)(array);
            for (i=0;i<num;i++)
                uchar_array[i] -= offset;
        break;
        case 2: ;
            unsigned short* short_array = (unsigned short*)(array);
            for (i=0;i<num;i++)
                short_array[i] -= offset;
        break;
        case 4: ;
            unsigned int* int_array = (unsigned int*)(array);
            for (i=0;i<num;i++) 
                int_array[i] -= offset;
        break;
        case 8: ; 
            unsigned long long* long_array = (unsigned long long*)(array);
            for (i=0;i<num;i++) 
                long_array[i] -= offset;           
        break;
        default:
            ffpmsg("ERROR: Non supported type size.");
        break;
    };
}

/* Pad the output array if the number of compressed values is less than the FITS array size
   @param array     The memory where the padded values should be written
   @param num       The number of values to write
   @param elem_size The size in bytes of a given value
   @param value     The value to use for padding
*/
void pad_with_value(char* array, 
                    int   num, 
                    int   elem_size, 
                    unsigned long long value)
{
    int i;
    switch (elem_size)
    {
        case 1: ;
            unsigned char* uchar_array = (unsigned char*)(array);
            for (i=0;i<num;i++)
                uchar_array[i] = value;
        break;
        case 2: ;
            unsigned short* short_array = (unsigned short*)(array);
            for (i=0;i<num;i++)
                short_array[i] = value;
           
        break;
        case 4: ;
            unsigned int* int_array = (unsigned int*)(array);
            for (i=0;i<num;i++) 
                int_array[i] = value;
                
        break;
        case 8: ; 
            unsigned long long* long_array = (unsigned long long*)(array);
            for (i=0;i<num;i++)
                long_array[i] = value;        
        break;
        default:
            printf("Non supported type size: %d\n", elem_size);
        break;
    };
}

/* Copy values from on array to another, changing their endian-ness along the way
   @param dest      The destination array
   @param src       The source array 
   @param num       The number of values to copy
   @param elem_size The size in bytes of a given value*/
void copy_swap(char*      dest, 
              const char* src, 
              int         num, 
              int         elem_size)
{
    int i;
    const char *pend = src + num*elem_size;
    for (const char *ptr = src; ptr<pend; ptr+=elem_size, dest+=elem_size)
        for (i=0;i<elem_size;i++)
            dest[i] = ptr[elem_size-i-1];
}

/* Add a new huffman encoding, a.k.a. symbol to the look-up-table
   @param lut   The top-level LUT node where to add the symbol
   @param sym   The decoded value
   @param nbits The number of bits occupied by the huffman encoding
   @param bits  The huffman-encoded bits */
void add_symbol(struct LUT**   lut, 
                unsigned short sym, 
                unsigned char  nbits, 
                size_t         bits)
{
    int i;
    // If no LUT available, create one
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
    
    /* From here, we deal with all the 8-bits leaves that contain the remaining code */
    int n_leaves = 1<<(8-nbits);
    for (i=0;i<n_leaves;i++) {
        const unsigned char key = bits | (i<<nbits);
        (*lut)[key].symbol = sym;
        (*lut)[key].is_leaf = 1;
        (*lut)[key].nbits   = nbits;
    }
}

/* Decompress a CTA-encoded buffer
   @param c           Input array
   @param input_len   Size of the input array in bytes
   @param array       Output array
   @param output_size Size of the output array in bytes
   @param col_type    Type of the FITS column in FITS enum value
   @param col_width   Width of the column in number of values
   @param col_zero    TZERO header key for that column
*/
int fits_ctadecomp(unsigned char* c,
                   unsigned long  input_len, 
                   unsigned char  array[], 
                   unsigned long  output_size, 
                   int            col_type,
                   int            col_width,
                   unsigned long long col_zero)
{
    unsigned long output_num = output_size / 2;
    unsigned char* max_input_address = c + input_len;

    struct BlockHeader* head = (struct BlockHeader*)(c);

    int i;

    c += sizeof(struct BlockHeader) + head->numProcs*sizeof(unsigned short);

    if (head->numProcs != 2)
    {
        ffpmsg("ERROR: Number of processings applied to CTA compression not 2 as expected.");
        return 1;
    }
    
    if (head->processings[0] != 0x3 || head->processings[1] != 0x2)
    {
        ffpmsg("ERROR: CTA compression processings not 0x3 and 0x2 as expected.");
        return 2;
    }

    /* Read the size of the compressed data */
    unsigned int buffer_length = 0;
    memcpy(&buffer_length, c, sizeof(buffer_length));
    c += sizeof(buffer_length);

    /* Read the number of 16 bits data this encoding represents. */
    unsigned long long data_count = 0;
    memcpy(&data_count, c, sizeof(data_count));
    c += sizeof(data_count);

    /* Read the number of symbols in the look-up table */
    unsigned long long sym_count = 0;
    memcpy(&sym_count, c, sizeof(sym_count));
    c += sizeof(sym_count);

    /* Read the look up table of the encoded symbols */
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
            ffpmsg("ERROR: Impossible Huffman encoding");
            return 1;
        }

        size_t bits = 0;
        memcpy(&bits, c, num_bytes);
        c += num_bytes;

        add_symbol(&lut, symbol, num_bits, bits);
    }

    /* if we have only one input symbol, fill up the output with it */
    if (sym_count == 1) {

        for (i=0;i<output_num;i++) {
            ((unsigned short*)(array))[i] = lut->symbol;
        }

        return 0;
    }
   
    /* take a temporary array to reshuffle the uncompressed data */
    unsigned short* tmp_array = malloc(data_count*sizeof(unsigned short));


    struct LUT top_lut;
    top_lut.lut = lut;
    top_lut.is_leaf = 0;
    struct LUT* clut = &top_lut; 
    unsigned char curbit = 0;
    unsigned long num_written = 0;

    /* Consume the compressed buffer, making sure we don't exhaust the input, nor write beyond the output buffer */
    while (c < max_input_address  && 
           num_written < data_count) {
        
        /* look at the current input byte. Take a two-bytes word as input to make sure that we can extract a full byte */
        const unsigned short *two = (unsigned short*)c;
        const unsigned char curbyte =  (*two >> curbit);
        
        /* if we end up nowhere, there is a problem with the encoding */
        if (!clut) {
            ffpmsg("Ended up with a NULL lut.");
            return 1;
        }

        /* otherwise, traverse the graph of symbols following the new input set of bits (curbyte) */
        clut = &(clut->lut[curbyte]);

        /* if we are not at a leaf yet, continue to traverse */
        if (!clut->is_leaf)
        {
            c++;
            continue;
        }

        /* if we've hit a leaf, assign the symbol to the output, and consume the corresponding bits */
        tmp_array[num_written++] = clut->symbol;

        curbit    += clut->nbits;
        /* and start over from the root of the graph (lut) */
        clut = &top_lut;
        if (curbit>=8){
            curbit %= 8;
            c++;
        }
    }

    /* Un-apply the subtraction of the previous element */
    for (i=1;i<data_count;i++)
        tmp_array[i] += tmp_array[i-1];

    int bytes_factor = 1;
    switch (col_type) {
        case TBIT:
            ffpmsg("TBITS are not handled by CTA compression...\n");
            return(-1);
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
        case TULONG:
        case TLONG:
        case TFLOAT:
            bytes_factor = 4;
        break;
        case TULONGLONG:
        case TLONGLONG:
        case TDOUBLE:
        case TCOMPLEX:
        case TDBLCOMPLEX:
            bytes_factor = 8;
        break;
        default:
            ffpmsg("ERROR: Got a default case in ctacomp.c.");
            return(-1);
        break;
    };

    /* Copy the temp array data to the output array */
    if (col_width == 1) {
        apply_offset((char*)(tmp_array), output_size / bytes_factor, bytes_factor, col_zero);
        copy_swap(array, (char*)(tmp_array), output_size / bytes_factor, bytes_factor);
    }
    else { /* Deal with arrays of possibly variable length*/
        int src_index = 0;
        int dst_index = 0;
        unsigned char* char_array = (unsigned char*)(tmp_array);
        while (src_index < (data_count*sizeof(unsigned short)-2)) { /* -2 because buffer size is rounded to 4 bytes sizes when writing*/ 
            
            int this_row_bytes = 0;
            memcpy(&this_row_bytes, &char_array[src_index], sizeof(this_row_bytes));
            src_index += sizeof(this_row_bytes);

            apply_offset( &(char_array[src_index]), this_row_bytes / bytes_factor, bytes_factor, col_zero);
            
            copy_swap(&(array[dst_index]), &(char_array[src_index]), this_row_bytes / bytes_factor, bytes_factor);
            dst_index += this_row_bytes;
            src_index += this_row_bytes;

            /* Pad any row shorter than the expected size with zeros */
            int missing_num = col_width - this_row_bytes/bytes_factor;
            pad_with_value(&(array[dst_index]), missing_num, bytes_factor, -col_zero);
            dst_index += missing_num*bytes_factor;
        }

        if (dst_index != output_size) {
            ffpmsg("WARNING WARNING: Written and expected sizes don't match");
            return(-1);
        }
    }

    free(tmp_array);

    erase_lut(lut);

    return 0;
}


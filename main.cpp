#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "algorithm.hpp"

int main()
{
    FILE* input = fopen ("input.txt", "r");
    if (input == NULL)
    {
        printf ("%s\n", strerror(errno));
        exit(-1);
    }
    const unsigned int block_size = 16; // in bytes, equals 128 bits
    /*
    char* buff = (char*) malloc (block_size * sizeof (char));
    if (buff == NULL)
    {
        printf ("malloc failed\n");
        exit(-1);
    }
    char symb = '0';
    for (unsigned int i = 0; i < block_size || symb != EOF; i++)
    {
        symb = fgetc(input);
        buff[i] = symb; 
    }
    */
    uint8_t buff [16] = {0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x00, 0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88};
    uint8_t key[32] =   {0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff, 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                         0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
    encrypt (buff, key, block_size);
    return 0;
}

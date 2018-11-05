#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

void encrypt (uint8_t* block, uint8_t* key, const unsigned int block_size);

/*
class Kuz
{
private:
    const unsigned int block_size = 16; //128 bits
    uint8_t C[32][block_size];
    uint8_t keys[10][block_size];
    uint8_t block[block_size];
public:
    Kuz ();
    void encrypt ();
    void decrypt ();
    ~Kuz ();
}
*/
#endif

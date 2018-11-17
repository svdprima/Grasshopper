#include <fstream>
#include <iostream>
#include "encryptor.hpp"

using std::fstream;
using std::cerr;
using std::cout;

void Encryptor::ReadText(const string& filename, bool as_encrypted) {
    Clear();
    fstream f(filename, std::ios_base::in);
    if (!f.is_open()) {
        cerr << "Cannot open file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    f.seekg(0, std::ios::end);
    text_size.u64 = f.tellg();
    f.seekg(0, std::ios::beg);

    uint8_t reserved_for_fsize = as_encrypted ? 0 : constant_size;

    uint64_t size = text_size.u64 + reserved_for_fsize;
    uint64_t num_blocks = size / Grasshopper::block_size + !!(size % Grasshopper::block_size);

    cout << "Number of blocks: " << num_blocks << "\n";

    text.resize(num_blocks);

    std::copy(text_size.u8, text_size.u8 + reserved_for_fsize, text.data()->begin());
    uint8_t* begin = text.data()->begin() + reserved_for_fsize;
    f.read(reinterpret_cast<char *>(begin), text_size.u64);
    if (!f) {
        cerr << "Cannot read all file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
    f.close();
}

void Encryptor::SaveText(const string& filename, bool as_decrypted) {
    fstream f(filename, std::ios_base::out);
    if (!f.is_open()) {
        cerr << "Cannot open file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    if (as_decrypted)
        std::copy(text.data()->begin(), text.data()->begin() + constant_size, text_size.u8);
    else
        text_size.u64 = text.size() * Grasshopper::block_size;

    uint64_t reserved_for_fsize = as_decrypted ? constant_size : 0;
    const uint8_t *begin = text.data()->begin() + reserved_for_fsize;
    f.write(reinterpret_cast<const char*>(begin), text_size.u64);
    if (!f) {
        cerr << "Cannot read all file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
    f.close();
}

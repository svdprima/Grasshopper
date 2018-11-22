#ifndef ENCRYPTOR_H
#define ENCRYPTOR_H

#include <string>
#include <vector>
#include "grasshopper.hpp"
#include "tables.hpp"

union UintConverter{
    uint64_t u64;
    uint8_t u8[8];
};


class Encryptor {
public:
    Encryptor() : text_size({0}), constant_size(sizeof(uint64_t)), mode(ECB) {}
    Encryptor(unsigned md) : text_size({0}), constant_size(sizeof(uint64_t)){mode = md;}
    void Encrypt(unsigned mode) {
      EncryptorImpl.Encrypt(text, InitialKey, mode);
    }
    void Decrypt(unsigned mode) {
      EncryptorImpl.Decrypt(text, InitialKey, mode);
    }
    void ReadText(const std::string& filename, bool as_encrypted = false);
    void SaveText(const std::string& filename, bool as_decrypted = false);
    void Clear() {
      text.clear();
      text_size.u64 = 0;
    }

private:
    Grasshopper EncryptorImpl;
    Grasshopper::Data text;
    UintConverter text_size;
    uint8_t constant_size;
    unsigned mode;

    Grasshopper::Key InitialKey {{0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
                                  0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                                  0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10,
                                  0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef}};
};

#endif

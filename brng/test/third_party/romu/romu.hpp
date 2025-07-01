#include <cstdint>

void romuQuad_init(uint64_t seed1, uint64_t seed2, uint64_t seed3,
                   uint64_t seed4);
uint64_t romuQuad_random();

void romuTrio_init(uint64_t seed1, uint64_t seed2, uint64_t seed3);
uint64_t romuTrio_random();

void romuDuo_init(uint64_t seed1, uint64_t seed2);
uint64_t romuDuo_random();

void romuDuoJr_init(uint64_t seed1, uint64_t seed2);
uint64_t romuDuoJr_random();

void romuQuad32_init(uint64_t seed1, uint64_t seed2, uint64_t seed3,
                     uint64_t seed4);
uint32_t romuQuad32_random();

void romuTrio32_init(uint64_t seed1, uint64_t seed2, uint64_t seed3);
uint32_t romuTrio32_random();

void romuMono32_init(uint32_t seed);
uint16_t romuMono32_random();

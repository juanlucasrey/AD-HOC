// Romu Pseudorandom Number Generators
//
// Copyright 2020 Mark A. Overton
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// ------------------------------------------------------------------------------------------------
//
// Website: romu-random.org
// Paper:   http://arxiv.org/abs/2002.11331
//
// Copy and paste the generator you want from those below.
// To compile, you will need to #include <stdint.h> and use the ROTL definition
// below.

#include <cstdint>

#define ROTL(d, lrot) ((d << (lrot)) | (d >> (8 * sizeof(d) - (lrot))))

//===== RomuQuad
//==================================================================================
//
// More robust than anyone could need, but uses more registers than RomuTrio.
// Est. capacity >= 2^90 bytes. Register pressure = 8 (high). State size = 256
// bits.

uint64_t wStateQuad, xStateQuad, yStateQuad, zStateQuad; // set to nonzero seed

void romuQuad_init(uint64_t seed1, uint64_t seed2, uint64_t seed3,
                   uint64_t seed4) {
    wStateQuad = seed1;
    xStateQuad = seed2;
    yStateQuad = seed3;
    zStateQuad = seed4;
}

uint64_t romuQuad_random() {
    uint64_t wp = wStateQuad, xp = xStateQuad, yp = yStateQuad, zp = zStateQuad;
    wStateQuad = 15241094284759029579u * zp; // a-mult
    xStateQuad = zp + ROTL(wp, 52);          // b-rotl, c-add
    yStateQuad = yp - xp;                    // d-sub
    zStateQuad = yp + wp;                    // e-add
    zStateQuad = ROTL(zStateQuad, 19);       // f-rotl
    return xp;
}

//===== RomuTrio
//==================================================================================
//
// Great for general purpose work, including huge jobs.
// Est. capacity = 2^75 bytes. Register pressure = 6. State size = 192 bits.

uint64_t xStateTrio, yStateTrio, zStateTrio; // set to nonzero seed

void romuTrio_init(uint64_t seed1, uint64_t seed2, uint64_t seed3) {
    xStateTrio = seed1;
    yStateTrio = seed2;
    zStateTrio = seed3;
}

uint64_t romuTrio_random() {
    uint64_t xp = xStateTrio, yp = yStateTrio, zp = zStateTrio;
    xStateTrio = 15241094284759029579u * zp;
    yStateTrio = yp - xp;
    yStateTrio = ROTL(yStateTrio, 12);
    zStateTrio = zp - yp;
    zStateTrio = ROTL(zStateTrio, 44);
    return xp;
}

//===== RomuDuo
//==================================================================================
//
// Might be faster than RomuTrio due to using fewer registers, but might
// struggle with massive jobs. Est. capacity = 2^61 bytes. Register pressure
// = 5. State size = 128 bits.

uint64_t xStateDuo, yStateDuo; // set to nonzero seed

void romuDuo_init(uint64_t seed1, uint64_t seed2) {
    xStateDuo = seed1;
    yStateDuo = seed2;
}

uint64_t romuDuo_random() {
    uint64_t xp = xStateDuo;
    xStateDuo = 15241094284759029579u * yStateDuo;
    yStateDuo = ROTL(yStateDuo, 36) + ROTL(yStateDuo, 15) - xp;
    return xp;
}

//===== RomuDuoJr
//================================================================================
//
// The fastest generator using 64-bit arith., but not suited for huge jobs.
// Est. capacity = 2^51 bytes. Register pressure = 4. State size = 128 bits.

uint64_t xStateDuoJr, yStateDuoJr; // set to nonzero seed

void romuDuoJr_init(uint64_t seed1, uint64_t seed2) {
    xStateDuoJr = seed1;
    yStateDuoJr = seed2;
}

uint64_t romuDuoJr_random() {
    uint64_t xp = xStateDuoJr;
    xStateDuoJr = 15241094284759029579u * yStateDuoJr;
    yStateDuoJr = yStateDuoJr - xp;
    yStateDuoJr = ROTL(yStateDuoJr, 27);
    return xp;
}

//===== RomuQuad32
//================================================================================
//
// 32-bit arithmetic: Good for general purpose use.
// Est. capacity >= 2^62 bytes. Register pressure = 7. State size = 128 bits.

uint32_t wStateQuad32, xStateQuad32, yStateQuad32,
    zStateQuad32; // set to nonzero seed

void romuQuad32_init(uint64_t seed1, uint64_t seed2, uint64_t seed3,
                     uint64_t seed4) {
    wStateQuad32 = seed1;
    xStateQuad32 = seed2;
    yStateQuad32 = seed3;
    zStateQuad32 = seed4;
}

uint32_t romuQuad32_random() {
    uint32_t wp = wStateQuad32, xp = xStateQuad32, yp = yStateQuad32,
             zp = zStateQuad32;
    wStateQuad32 = 3323815723u * zp;      // a-mult
    xStateQuad32 = zp + ROTL(wp, 26);     // b-rotl, c-add
    yStateQuad32 = yp - xp;               // d-sub
    zStateQuad32 = yp + wp;               // e-add
    zStateQuad32 = ROTL(zStateQuad32, 9); // f-rotl
    return xp;
}

//===== RomuTrio32
//===============================================================================
//
// 32-bit arithmetic: Good for general purpose use, except for huge jobs.
// Est. capacity >= 2^53 bytes. Register pressure = 5. State size = 96 bits.

uint32_t xStateTrio32, yStateTrio32, zStateTrio32; // set to nonzero seed

void romuTrio32_init(uint64_t seed1, uint64_t seed2, uint64_t seed3) {
    xStateTrio32 = seed1;
    yStateTrio32 = seed2;
    zStateTrio32 = seed3;
}

uint32_t romuTrio32_random() {
    uint32_t xp = xStateTrio32, yp = yStateTrio32, zp = zStateTrio32;
    xStateTrio32 = 3323815723u * zp;
    yStateTrio32 = yp - xp;
    yStateTrio32 = ROTL(yStateTrio32, 6);
    zStateTrio32 = zp - yp;
    zStateTrio32 = ROTL(zStateTrio32, 22);
    return xp;
}

//===== RomuMono32
//===============================================================================
//
// 32-bit arithmetic: Suitable only up to 2^26 output-values. Outputs 16-bit
// numbers. Fixed period of (2^32)-47. Must be seeded using the romuMono32_init
// function. Capacity = 2^27 bytes. Register pressure = 2. State size = 32 bits.

uint32_t stateMono32;

void romuMono32_init(uint32_t seed) {
    stateMono32 = (seed & 0x1fffffffu) + 1156979152u; // Accepts 29 seed-bits.
}

uint16_t romuMono32_random() {
    uint16_t result = stateMono32 >> 16;
    stateMono32 *= 3611795771u;
    stateMono32 = ROTL(stateMono32, 12);
    return result;
}

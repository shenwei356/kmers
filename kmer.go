// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//b
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package kmers

import (
	"bytes"
	"errors"
	"math/bits"
)

// ErrIllegalBase means that base beyond IUPAC symbols are  detected.
var ErrIllegalBase = errors.New("kmers: illegal base")

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("kmers: k-mer size (1-32) overflow")

// ErrPositionOverflow means i >= K
var ErrPositionOverflow = errors.New("kmers: base position (0-based) overflow")

// ErrLengthOverflow means the length n > K
var ErrLengthOverflow = errors.New("kmers: base position (0-based) overflow")

// ErrCodeOverflow means the encode interger is bigger than 4^k.
var ErrCodeOverflow = errors.New("kmers: code value overflow")

// ErrKMismatch means K size mismatch.
var ErrKMismatch = errors.New("kmers: K mismatch")

// slice is much faster than switch and map.
var base2bit = [256]uint64{
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 1, 1, 0, 4, 4, 2, 0, 4, 4, 2, 4, 0, 0, 4,
	4, 4, 0, 1, 3, 3, 0, 0, 4, 1, 4, 4, 4, 4, 4, 4,
	4, 0, 1, 1, 0, 4, 4, 2, 0, 4, 4, 2, 4, 0, 0, 4,
	4, 4, 0, 1, 3, 3, 0, 0, 4, 1, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
}

// var base2bit []uint64

// MaxCode is the maxinum interger for all Ks.
var MaxCode []uint64

func init() {
	MaxCode = make([]uint64, 33)
	for i := 1; i <= 32; i++ {
		MaxCode[i] = (1 << uint(1<<i)) - 1 // (1<<(i*2)) - 1
	}
}

// Encode converts byte slice to bits.
//
// Codes:
//
//	A    0b00
//	C    0b01
//	G    0b10
//	T    0b11
//
// For degenerate bases, only the first base is kept.
//
//	M       AC     A
//	V       ACG    A
//	H       ACT    A
//	R       AG     A
//	D       AGT    A
//	W       AT     A
//	S       CG     C
//	B       CGT    C
//	Y       CT     C
//	K       GT     G
//	N       ACGT   A
func Encode(kmer []byte) (code uint64, err error) {
	if len(kmer) == 0 || len(kmer) > 32 {
		return 0, ErrKOverflow
	}

	var v uint64
	for _, b := range kmer {
		code <<= 2
		v = base2bit[b]
		// if v > 3 {
		if v == 4 {
			return code, ErrIllegalBase
		}
		code |= v
	}
	return code, nil
}

// ErrNotConsecutiveKmers means the two k-mers are not consecutive.
var ErrNotConsecutiveKmers = errors.New("kmers: not consecutive k-mers")

// MustEncodeFromFormerKmer encodes from former the k-mer,
// assuming the k-mer and leftKmer are both OK.
func MustEncodeFromFormerKmer(kmer []byte, leftKmer []byte, leftCode uint64) (uint64, error) {
	v := base2bit[kmer[len(kmer)-1]]
	// if v > 3 {
	if v == 4 {
		return leftCode, ErrIllegalBase
	}
	// retrieve lower (k-1)*2 bits and << 2, and then add v
	return (leftCode&((1<<(uint(len(kmer)-1)<<1))-1))<<2 | v, nil
}

// EncodeFromFormerKmer encodes from the former k-mer, inspired by ntHash
func EncodeFromFormerKmer(kmer []byte, leftKmer []byte, leftCode uint64) (uint64, error) {
	if len(kmer) == 0 {
		return 0, ErrKOverflow
	}
	if len(kmer) != len(leftKmer) {
		return 0, ErrKMismatch
	}
	if !bytes.Equal(kmer[0:len(kmer)-1], leftKmer[1:]) {
		return 0, ErrNotConsecutiveKmers
	}
	return MustEncodeFromFormerKmer(kmer, leftKmer, leftCode)
}

// MustEncodeFromLatterKmer encodes from the latter k-mer,
// assuming the k-mer and rightKmer are both OK.
func MustEncodeFromLatterKmer(kmer []byte, rightKmer []byte, rightCode uint64) (uint64, error) {
	v := base2bit[kmer[0]]
	// if v > 3 {
	if v == 4 {
		return rightCode, ErrIllegalBase
	}

	return v<<(uint(len(kmer)-1)<<1) | rightCode>>2, nil
}

// EncodeFromLatterKmer encodes from the former k-mer.
func EncodeFromLatterKmer(kmer []byte, rightKmer []byte, rightCode uint64) (uint64, error) {
	if len(kmer) == 0 {
		return 0, ErrKOverflow
	}
	if len(kmer) != len(rightKmer) {
		return 0, ErrKMismatch
	}
	if !bytes.Equal(rightKmer[0:len(kmer)-1], kmer[1:len(rightKmer)]) {
		return 0, ErrNotConsecutiveKmers
	}
	return MustEncodeFromLatterKmer(kmer, rightKmer, rightCode)
}

// Reverse returns code of the reversed sequence.
func Reverse(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code & 3)
	// 	code >>= 2
	// }
	// return

	// https: //www.biostars.org/p/113640, with a little modification
	c = code
	c = (c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2
	c = (c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4
	c = (c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8
	c = (c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16
	c = (c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32
	return c >> ((32 - k) << 1)
}

// MustReverse is similar to Reverse, but does not check k.
func MustReverse(code uint64, k int) (c uint64) {
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code & 3)
	// 	code >>= 2
	// }
	// return

	// https: //www.biostars.org/p/113640, with a little modification
	c = code
	c = (c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2
	c = (c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4
	c = (c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8
	c = (c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16
	c = (c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32
	return c >> ((32 - k) << 1)
}

// Complement returns code of complement sequence.
func Complement(code uint64, k int) uint64 {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	return code ^ ((1 << uint(k<<1)) - 1)
}

// MustComplement is similar to Complement, but does not check k.
func MustComplement(code uint64, k int) uint64 {
	return code ^ ((1 << uint(k<<1)) - 1)
}

// RevComp returns code of reverse complement sequence.
func RevComp(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code&3 ^ 3)
	// 	code >>= 2
	// }
	// return

	// https://www.biostars.org/p/113640/#9474334
	c = ^code
	c = (c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2
	c = (c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4
	c = (c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8
	c = (c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16
	c = (c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32
	return c >> ((32 - k) << 1)
}

// MustRevComp is similar to RevComp, but does not check k.
func MustRevComp(code uint64, k int) (c uint64) {
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code&3 ^ 3)
	// 	code >>= 2
	// }
	// return

	// https://www.biostars.org/p/113640/#9474334
	c = ^code
	c = (c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2
	c = (c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4
	c = (c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8
	c = (c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16
	c = (c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32
	return c >> ((32 - k) << 1)
}

// Canonical returns code of its canonical kmer.
func Canonical(code uint64, k int) uint64 {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}

	var rc uint64
	// c := code
	// for i := 0; i < k; i++ {
	// 	rc = (rc << 2) | (c&3 ^ 3)
	// 	c >>= 2
	// }

	// https://www.biostars.org/p/113640/#9474334
	c := ^code
	c = (c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2
	c = (c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4
	c = (c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8
	c = (c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16
	c = (c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32
	rc = c >> ((32 - k) << 1)

	if rc < code {
		return rc
	}
	return code
}

// MustCanonical is similar to Canonical, but does not check k.
func MustCanonical(code uint64, k int) uint64 {
	var rc uint64
	// c := code
	// for i := 0; i < k; i++ {
	// 	rc = (rc << 2) | (c&3 ^ 3)
	// 	c >>= 2
	// }

	// https://www.biostars.org/p/113640/#9474334
	c := ^code
	c = (c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2
	c = (c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4
	c = (c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8
	c = (c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16
	c = (c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32
	rc = c >> ((32 - k) << 1)

	if rc < code {
		return rc
	}
	return code
}

// just use: (1413956417 >> (x << 3)) & 255
//
// chr T        G        C        A
// 0b 01010100 01000111 01000011 01000001  1413956417

// bit2base is for mapping bit to base.
var bit2base = [4]byte{'A', 'C', 'G', 'T'}

// bit2str is for output bits string
var bit2str = [4]string{"00", "01", "10", "11"}

// Decode converts the code to original seq
func Decode(code uint64, k int) []byte {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	if code > MaxCode[k] {
		panic(ErrCodeOverflow)
	}
	kmer := make([]byte, k)
	for i := 0; i < k; i++ {
		kmer[k-1-i] = bit2base[code&3]

		// it's slower than bit2base[code&3], according to the test.
		// kmer[k-1-i] = byte(uint64(1413956417>>((code&3)<<3)) & 255)

		code >>= 2
	}
	return kmer
}

// MustDecode is similar to Decode, but does not check k and code.
func MustDecode(code uint64, k int) []byte {
	kmer := make([]byte, k)
	for i := 0; i < k; i++ {
		kmer[k-1-i] = bit2base[code&3]

		// it's slower than bit2base[code&3], according to the test.
		// kmer[k-1-i] = byte(uint64(1413956417>>((code&3)<<3)) & 255)

		code >>= 2
	}
	return kmer
}

// BaseAt returns the base in pos i (0-based).
func BaseAt(code uint64, k int, i int) uint8 {
	if i < 0 || i >= k {
		panic(ErrPositionOverflow)
	}
	return uint8(code >> ((k - i - 1) << 1) & 3)
}

// MustBaseAt returns the base in pos i (0-based).
func MustBaseAt(code uint64, k int, i int) uint8 {
	return uint8(code >> ((k - i - 1) << 1) & 3)
}

// Prefix returns the first n bases. n needs to be > 0.
// The length of the prefix is n.
func Prefix(code uint64, k int, n int) uint64 {
	if n < 1 || n > k {
		panic(ErrLengthOverflow)
	}
	return code >> ((k - n) << 1)
}

// MustPrefix returns the first n bases. n needs to be > 0.
// The length of the prefix is n.
func MustPrefix(code uint64, k int, n int) uint64 {
	return code >> ((k - n) << 1)
}

// Suffix returns the suffix starting from position i (0-based).
// The length of the suffix is k - i.
func Suffix(code uint64, k int, i int) uint64 {
	if i < 0 || i >= k {
		panic(ErrPositionOverflow)
	}
	return code & (1<<((k-i)<<1) - 1)
}

// MustSuffix returns the suffix starting from position i (0-based).
// The length of the suffix is k - i.
func MustSuffix(code uint64, k int, i int) uint64 {
	return code & (1<<((k-i)<<1) - 1)
}

// LongestPrefix returns the length of the longest prefix.
func LongestPrefix(code1, code2 uint64, k1, k2 int) int {
	if k1 <= 0 || k1 > 32 || k2 <= 0 || k2 > 32 {
		panic(ErrKOverflow)
	}

	var d int
	if k1 >= k2 { // most of the cases
		code1 = code1 >> ((k1 - k2) << 1)
		d = 32 - k2
	} else {
		code2 = code2 >> ((k2 - k1) << 1)
		d = 32 - k1
	}
	return bits.LeadingZeros64(code1^code2)>>1 - d
}

// MustLongestPrefix returns the length of the longest prefix.
func MustLongestPrefix(code1, code2 uint64, k1, k2 int) int {
	var d int
	if k1 >= k2 { // most of the cases
		code1 = code1 >> ((k1 - k2) << 1)
		d = 32 - k2
	} else {
		code2 = code2 >> ((k2 - k1) << 1)
		d = 32 - k1
	}
	return bits.LeadingZeros64(code1^code2)>>1 - d
}

// HasPrefix check if a k-mer has a prefix
func HasPrefix(code uint64, prefix uint64, k1, k2 int) bool {
	if k1 <= 0 || k1 > 32 || k2 <= 0 || k2 > 32 {
		panic(ErrKOverflow)
	}

	if k1 < k2 {
		return false
	}
	return code>>((k1-k2)<<1) == prefix
}

// MustHasPrefix check if a k-mer has a prefix
func MustHasPrefix(code uint64, prefix uint64, k1, k2 int) bool {
	if k1 < k2 {
		return false
	}
	return code>>((k1-k2)<<1) == prefix
}

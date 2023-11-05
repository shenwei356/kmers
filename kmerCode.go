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
)

// KmerCode is a struct representing a k-mer in 64-bits.
type KmerCode struct {
	Code uint64
	K    int
}

// NewKmerCode returns a new KmerCode struct from byte slice.
func NewKmerCode(kmer []byte) (KmerCode, error) {
	code, err := Encode(kmer)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// NewKmerCodeFromFormerOne computes KmerCode from the Former consecutive k-mer.
func NewKmerCodeFromFormerOne(kmer []byte, leftKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := EncodeFromFormerKmer(kmer, leftKmer, preKcode.Code)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// NewKmerCodeMustFromFormerOne computes KmerCode from the Former consecutive k-mer,
// assuming the k-mer and leftKmer are both OK.
func NewKmerCodeMustFromFormerOne(kmer []byte, leftKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := MustEncodeFromFormerKmer(kmer, leftKmer, preKcode.Code)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// Equal checks wether two KmerCodes are the same.
func (kcode KmerCode) Equal(kcode2 KmerCode) bool {
	return kcode.K == kcode2.K && kcode.Code == kcode2.Code
}

// Rev returns KmerCode of the reverse sequence.
func (kcode KmerCode) Rev() KmerCode {
	return KmerCode{MustReverse(kcode.Code, kcode.K), kcode.K}
}

// Comp returns KmerCode of the complement sequence.
func (kcode KmerCode) Comp() KmerCode {
	return KmerCode{MustComplement(kcode.Code, kcode.K), kcode.K}
}

// RevComp returns KmerCode of the reverse complement sequence.
func (kcode KmerCode) RevComp() KmerCode {
	return KmerCode{MustRevComp(kcode.Code, kcode.K), kcode.K}
}

// Canonical returns its canonical kmer
func (kcode KmerCode) Canonical() KmerCode {
	rcKcode := kcode.RevComp()
	if rcKcode.Code < kcode.Code {
		return rcKcode
	}
	return kcode
}

// Bytes returns k-mer in []byte.
func (kcode KmerCode) Bytes() []byte {
	return Decode(kcode.Code, kcode.K)
}

// String returns k-mer in string
func (kcode KmerCode) String() string {
	return string(Decode(kcode.Code, kcode.K))
}

// BitsString returns code to string
func (kcode KmerCode) BitsString() string {
	var buf bytes.Buffer
	for _, b := range Decode(kcode.Code, kcode.K) {
		buf.WriteString(bit2str[base2bit[b]])
	}
	return buf.String()
}

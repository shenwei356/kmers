# kmers

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/kmers.svg)](https://pkg.go.dev/github.com/shenwei356/kmers)


This package provides manipulations for bit-packed k-mers (k<=32, encoded in `uint64`).

Related projects:

- [unik](https://github.com/shenwei356/unik) provides k-mer serialization methods for this package.
- [unikmer](https://github.com/shenwei356/unikmer), a toolkit for nucleic acid k-mer analysis,
 including set operations on k-mers optional with TaxIDs.
- [sketches](https://pkg.go.dev/github.com/shenwei356/bio/sketches) provides generators/iterators for k-mer sketches 
([Minimizer](https://academic.oup.com/bioinformatics/article/20/18/3363/202143),
 [Scaled MinHash](https://f1000research.com/articles/8-1006),
 [Closed Syncmers](https://peerj.com/articles/10805/)).

## Benchmark

CPU: AMD Ryzen 7 2700X Eight-Core Processor, 3.7 GHz

    $ go test . -bench=Bench* -benchmem \
        | grep Bench \
        | perl -pe 's/\s\s+/\t/g' \
        | csvtk cut -Ht -f 1,3-5 \
        | csvtk add-header -t -n test,time,memory,allocs \
        | csvtk pretty -t -r
 
                                          test           time     memory        allocs
    ------------------------------------------   ------------   --------   -----------
                         BenchmarkEncodeK32-16    19.67 ns/op     0 B/op   0 allocs/op
           BenchmarkEncodeFromFormerKmerK32-16    7.692 ns/op     0 B/op   0 allocs/op
       BenchmarkMustEncodeFromFormerKmerK32-16    2.008 ns/op     0 B/op   0 allocs/op
                         BenchmarkDecodeK32-16    80.73 ns/op    32 B/op   1 allocs/op
                     BenchmarkMustDecodeK32-16    76.93 ns/op    32 B/op   1 allocs/op
    
                            BenchmarkRevK32-16    3.617 ns/op     0 B/op   0 allocs/op
                           BenchmarkCompK32-16   0.7999 ns/op     0 B/op   0 allocs/op
                        BenchmarkRevCompK32-16    3.814 ns/op     0 B/op   0 allocs/op
                       BenchmarkCannonalK32-16    4.147 ns/op     0 B/op   0 allocs/op

## Support

Please [open an issue](https://github.com/shenwei356/kmers/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/kmers/blob/master/LICENSE)

## History

This package was originally maintained in [unikmer](https://github.com/shenwei356/unikmer).

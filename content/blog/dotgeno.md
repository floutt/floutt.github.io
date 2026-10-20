Title: Unpacking a Tricky Format
Date: 2026-10-20
Author: Oluwatosin Olayinka
Slug: dotgeno

I've always been interested in human diversity, whether this comes in the form of differing societal norms, subsistence strategies, languages, or food. One aspect of this diversity that has always intrigued me is human genetic diversity - both within and between human populations. While we humans are all close genetic relatives of each other, the differences that we do have highlight the unique stories and challenges of our particular ancestors. For example, while most humans lose the ability to digest lactose after weaning, Northern Europeans and certain populations in East Africa and the Sahel have high rates of alleles which allow them to digest lactose for the entirety of their life. The high prevalence of these alleles is likely an adaptation to a diet uniquely high in unfermented dairy products.

Over the past two decades, developments in molecular biology have revolutionized our ability to study human genetic diversity. Not only can we now get a genetic picture of contemporary individuals across their entire genomes, we can also do so for people in the past. Thanks to the concurrent development of high-throughput genetic assays (whole-genome sequencing, DNA microarrays, etc.)<sup>1,2</sup> and more sophisticated tools to extract DNA from ancient human (and animal) fossils<sup>3</sup>, we can now read the genomes of _thousands_ of modern AND ancient individuals from all over the world. This has been revolutionary for the field of anthropology, where old hypotheses can now be tested at the molecular level. For example, many have speculated on archeological<sup>4,5</sup> and linguistic<sup>6,7</sup> bases of pre-Columbian contact between ancient Polynesian and South American populations. Thanks to the genomes of modern<sup>8</sup> AND ancient<sup>9</sup> Rapa Nui people we have a definitive answer - YES. Both modern and ancient Rapa Nui individuals show signals of admixture with native South Americans dating to hundreds of years before European contact.

One of the best resources for studying the genetics of ancient and modern human populations is the Allen Ancient DNA Resource (AADR)<sup>10</sup>. This is a curated dataset of thousands of samples from hundreds of different studies which have been aggregated together into a single set of files. This has lowered the barrier of entry for hobbyists interested in the topic of human population genetics, thus democratizing one aspect of the science. Despite all of its benefits, usage of the AADR is hampered by its use of bespoke, sparsely-supported file formats. 

## Struggles with the PACKEDANCESTRYMAP and TGENO formats
Genomic data from the AADR is, depending on the version, distributed in either the `PACKEDANCESTRYMAP` or `TGENO` file format, with the latter being used in more recent versions. Both of these file types come with corresponding `.snp` and `.ind` files storing metadata about genetic variants and individuals respectively. The `PACKEDANCESTRYMAP` and `TGENO` formats are closely related but differ in minor, but important, ways (**Figure 1**). In `PACKEDANCESTRYMAP` files, data is read on a "SNP-by-SNP" basis, where each row of data contains allelic dosages for a single variant across all individuals in the file. The `TGENO` format is the transposed version of this where data is read on an "individual-by-indivdiual" basis, with each row storing allelic dosages for a single individual across all variants in the file. 

![]({static}/images/dotgeno/schematic_bin_files.png)

: **Figure 1**: Schematic overview of the `PACKEDANCESTRYMAP` and `TGENO` file formats showing high-level similarities and differences between the two.

The `PACKEDANCESTRYMAP` format is widely used in population genetics software. For instance, it is the preferred input format for `smartpca`, `AdmixTools`, `ADMIXTOOLS 2`, and in a wide variety of software for admixture dating such as `DATES`, `ALDER`, `MALDER`, and `ROLLOF`. Despite its widespread usage in this domain, tools for directly manipulating this file type are rare, if not non-existent. This is in stark contrast to other genomic file types such as the Variant Call (`.vcf`), PLINK1 (`.bed`), and PLINK2 (`.pgen`) formats which are all widely supported and have robust tools for input fitering. This limitation becomes problematic when performing any analysis that requires complex filtering of the data. Want to perform a test of sex bias using alleles from the X-chromosome? Want to identify the allele frequencies of a particular set of alleles? Want to filter variants by missingness? Want to merge the AADR with another dataset? All of these tasks are non-trivial to perform with `PACKEDANCESTRYMAP` and `TGENO` files, so much so that the standard advice is to first convert the `PACKEDANCESTRYMAP` or `TGENO` file to a PLINK1 file, perform the necessary tasks with PLINK1, and then convert it BACK to its original format (**Figure 2**). This process of going back and forth between different intermediate file types is a waste of compute, memory, and storage space - precious resources for the hobbyist who lacks access to high-performance computing clusters.

![]({static}/images/dotgeno/genarchivist_screenshot.png)

: **Figure 2**: On the popular population genetics forum [GenArchivist](https://genarchivist.net/index.php), the general advice for merging external datasets with the AADR is to first [convert the `PACKEDANCESTRYMAP` file to a `PLINK1` file](https://genarchivist.net/showthread.php?tid=1279&page=3).

## Intoducing dotgenolib
In response to the frustrastions that come with these file formats, I took matters into my own hands. I developed [`dotgenolib`](https://github.com/floutt/dotgenolib/), a C library for reading and writing `PACKEDANCESTRYMAP`, `EIGENSTRAT`, and `TGENO` files. Not only is this library fast, it also comes with an easy to use, well-documented API and requires no external dependencies in a standard POSIX environment. I hope that by releasing this library to the public, hobbyists and researchers alike will have a much easier time working with this file type (and thus the AADR) than I have had in the past.

In this article I will go over the basics of how the library works and provide some real-world applications of the library to tasks and challeneges I've had. In cases where existing methods were available, comparisons were made between the `dotgenolib`-based implementation and the already existing one. In all cases, the `dotgenolib`-based implementation showed superior performance metrics to its standard alternative.

> **A note on compilation**: all of the code snippets in this post can be compiled with `gcc ${source_file} -o ${outbin} -ldotgeno -lm -Wl,-rpath,/usr/local/lib` where `${source_file}` is the path to the file with the code snippet in it and `${outbinm}` is the path to the executable. 

### How it works
`dotgenolib` is a simple library designed to be an easy to use, low-level library for reading and writing `PACKEDANCESTRYMAP`, `EIGENSTRAT`, and `TGENO` files. The library provides a robust set of functions for dealing with `.ind` and `.snp` metadata files and manipulating binary and text-based genotype files. The basics of the library can be illustrated by using it to convert a `PACKEDANCESTRYMAP` file to a simple SNP×sample table:

```c
#include <stdio.h>
#include <stdlib.h>
#include <dotgeno.h>

int main(int argc, char* argv[]) {
    if(argc != 4) {
        fprintf(stderr, "ERROR: not enough arguments.\n");
        return EXIT_FAILURE;
    }

    char* snp_file = argv[1];
    char* ind_file = argv[2];
    char* geno_file = argv[3];

    snp_data snps = read_snp_file(snp_file);
    ind_data inds = read_ind_file(ind_file);
    pam_file_reader pfr = pam_file_reader_init(geno_file, &snps, &inds);

    /* write output header */
    for(size_t i = 0; i < inds.length; i++) {
            printf("\t%s", inds.ind_id[i]);
    }
    printf("\n");

    /* write SNP records */ 
    /* for PACKEDANCESTRYMAP and TGENO files, metadata header must be read
    * before reading genotype records */
    read_pam_header(&pfr);
    uint8_t* dosages = NULL; // array for storing genotype dosages
    size_t snp_idx = 0;
    while(dosages = read_pam_record(&pfr)) {        
            printf("%s", snps.var_id[snp_idx]);
            if(dosages == NULL) { break; }
            for(size_t i = 0; i < inds.length; i++) {
                    printf("\t%hhu", dosages[i]);
            }
            printf("\n");
            free(dosages);
            snp_idx++;
    }
    /* free data and close file streams */
    free_snp_data(&snps);
    free_ind_data(&inds);
    close_pam_file_reader(&pfr);
}
```

When we run the above program with example [`.snp`](https://github.com/DReichLab/AdmixTools/blob/master/convertf/example.snp), [`.ind`](https://github.com/DReichLab/AdmixTools/blob/master/convertf/example.ind), and [`.packedancestrymapgeno`](https://github.com/DReichLab/AdmixTools/blob/master/convertf/example.packedancestrymapgeno) files we get the following output: 

```
        SAMPLE0 SAMPLE1 SAMPLE2 SAMPLE3 SAMPLE4
rs0000  1       1       1       0       0
rs1111  0       1       2       1       2
rs2222  2       1       1       0       1
rs3333  0       0       1       2       2
rs4444  2       1       1       0       0
rs5555  0       0       1       1       1
rs6666  2       2       1       1       0
```

The library also contains functions for filtering and writing genotype files. We can use them to, for example, filter a `PACKEDANCESTRYMAP` file to only female samples:

```c
#include <stdio.h>
#include <stdlib.h>
#include <sys/queue.h> // dotgenolib uses TAILQ linked lists
#include <dotgeno.h>

int main(int argc, char* argv[]) {
	if(argc != 7) {
		fprintf(stderr, "ERROR: not enough arguments.\n");
		return EXIT_FAILURE;
	}

	char* in_snp_file = argv[1];
	char* in_ind_file = argv[2];
	char* in_geno_file = argv[3];
	char* out_snp_file = argv[4];
	char* out_ind_file = argv[5];
	char* out_geno_file = argv[6];

	snp_data in_snp = read_snp_file(in_snp_file);
	ind_data in_ind = read_ind_file(in_ind_file);

	/* linked list storing indices of female samples */
	struct idx_head females_list;
	TAILQ_INIT(&females_list);
	get_multiple_sex(&in_ind, "F", &females_list);

	/* write output snp and ind files */
	ind_data out_ind;
	filter_ind_data(&in_ind, &out_ind, &females_list);
	write_ind_data(&out_ind, out_ind_file);
	write_snp_data(&in_snp, out_snp_file); // SNPs not filtered

	/* reader and writer */
	pam_file_reader p_in = pam_file_reader_init(in_geno_file, &in_snp, &in_ind);
	pam_file_writer p_out = pam_file_writer_init(out_geno_file, &in_snp, &out_ind);
	
	// handle header
	read_pam_header(&p_in);
	write_pam_header(&p_out, &in_snp, &out_ind);

	uint8_t* in_record = NULL;
	uint8_t* out_record = (uint8_t*)malloc(sizeof(uint8_t) * out_ind.length);
	
	while(in_record = read_pam_record(&p_in)) {
		struct idx_node* np;
		size_t i = 0;
		TAILQ_FOREACH(np, &females_list, nodes) {
			out_record[i] = in_record[np->idx];
			i++;	
		}
		write_pam_record(&p_out, out_record);
		free(in_record);
	}
	
    free(out_record);
	free_idx_list(&females_list);
	free_ind_data(&in_ind);
	free_ind_data(&out_ind);
	free_snp_data(&in_snp);
	close_pam_file_reader(&p_in);
	close_pam_file_writer(&p_out);
}
```

When we run the above program on the same set of files and run the SNP×sample table code with the output `PACKEDANCESTRYMAP` files as input we get the following result:

```
        SAMPLE0 SAMPLE2 SAMPLE4
rs0000  1       1       0
rs1111  0       2       2
rs2222  2       1       1
rs3333  0       1       2
rs4444  2       1       0
rs5555  0       1       1
rs6666  2       1       0
```

## Real world use cases
These toy examples are fun and cool but how does the library perform in the wild? In what ways can and has it been used? In the next few sections, I will provide some example use cases of the library and, when possible, benchmark the `dotgenolib`-based implementations with pre-existing methods. 

> **A note on benchmarks**: all benchmarks were performed on a Thinkpad T14 Gen 6 laptop with 16 GB of RAM and an Intel® Core™ Ultra 5 225U.

### Converting between EIGENSTRAT and PACKEDANCESTRYMAP formats
`dotgenolib` can be easily be used to convert between the text-based `EIGENSTRAT` and the binary `PACKEDANCESTRYMAP` format, requiring very few lines of code to do so. Even with this simplicity, these conversions greatly outperform the standard `convertf` workkflow, with the `dotgenolib` implementations running up to 55× faster and with 35× less peak memory usage than `convertf` (**Table 1**). Converting between `PACKEDANCESTRYMAP` and `EIGENSTRAT` formats with this library is simple, as illustrated with the following code:

```c
#include <stdio.h>
#include <stdlib.h>
#include <dotgeno.h>

int main(int argc, char* argv[]) {
	if(argc != 7) {
		fprintf(stderr, "ERROR: not enough arguments.\n");
		return EXIT_FAILURE;
	}

	char* in_snp_file = argv[1];
	char* in_ind_file = argv[2];
	char* in_geno_file = argv[3];
	char* out_snp_file = argv[4];
	char* out_ind_file = argv[5];
	char* out_geno_file = argv[6];

	snp_data snp = read_snp_file(in_snp_file);
	ind_data ind = read_ind_file(in_ind_file);
	/* output ind and snp are the same */
	write_snp_data(&snp, out_snp_file);
	write_ind_data(&ind, out_ind_file);
	/* reader and writer objects */
	pam_file_reader pfr = pam_file_reader_init(in_geno_file, &snp, &ind);
	egn_file_writer efw = egn_file_writer_init(out_geno_file, &snp, &ind);
	read_pam_header(&pfr);
	uint8_t* dosages = NULL; // array for storing genotype dosages
	while(dosages = read_pam_record(&pfr)) {
		write_egn_record(&efw, dosages);	
		free(dosages);
	}
	/* free data and close file streams */
	free_snp_data(&snp);
	free_ind_data(&ind);
	close_pam_file_reader(&pfr);
	close_egn_file_writer(&efw);
}
```

In this library, the general process of reading and writing both file types is similar so converting the above code from `PACKEDANCESTRMAP`→`EIGENSTRAT` to `EIGENSTRAT`→`PACKEDANCESTRMAP` conversion requires minimal changes. Outside of some minor changes in header logic, swapping `pam` and `egn` with each other in the code gets most of the job done. The code for this conversion can be seen below

```c
#include <stdio.h>
#include <stdlib.h>
#include <dotgeno.h>

int main(int argc, char* argv[]) {
	if(argc != 7) {
		fprintf(stderr, "ERROR: not enough arguments.\n");
		return EXIT_FAILURE;
	}

	char* in_snp_file = argv[1];
	char* in_ind_file = argv[2];
	char* in_geno_file = argv[3];
	char* out_snp_file = argv[4];
	char* out_ind_file = argv[5];
	char* out_geno_file = argv[6];

	snp_data snp = read_snp_file(in_snp_file);
	ind_data ind = read_ind_file(in_ind_file);
	/* output ind and snp are the same */
	write_snp_data(&snp, out_snp_file);
	write_ind_data(&ind, out_ind_file);
	/* reader and writer objects */
	egn_file_reader efr = egn_file_reader_init(in_geno_file, &snp, &ind);
	pam_file_writer pfw = pam_file_writer_init(out_geno_file, &snp, &ind);
	write_pam_header(&pfw, &snp, &ind);
	uint8_t* dosages = NULL; // array for storing genotype dosages
	while(dosages = read_egn_record(&efr)) {
		write_pam_record(&pfw, dosages);	
		free(dosages);
	}
	/* free data and close file streams */
	free_snp_data(&snp);
	free_ind_data(&ind);
	close_egn_file_reader(&efr);
	close_pam_file_writer(&pfw);
}
```

For benchmarking, `PACKEDANCESTRYMAP` and `EIGENSTRAT` files with N SNPs and N individuals were generated for increasing values of N. From here, conversions between `PACKEDANCESTRYMAP` and `EIGENSTRAT` formats (and vice versa) were performed using both the aforementioned `dotgenolib`-based code and the standard `convertf` based method. For each conversion, peak memory usage (maximum RSS) and run time were compared between the `dotgenolib`-based code and `convertf`. For both conversions there was a very strong reduction in memory usage and run time in the `dotgenolib` based conversions (**Figure 3**, **Table 1**). 

![]({static}/images/dotgeno/mem_speed.png)

: **Figure 3**: Comparison of execution time and maximum resident set size (RSS) between `dotgenolib`-based programs and `convertf` for `EIGENSTRAT` to `PACKEDANCESTRYMAP` (EGN-TO-PAM) and `PACKEDANCESTRYMAP` to `EIGENSTRAT` (PAM-TO-EGN) conversions.


![]({static}/images/dotgeno/mem_speed_alt.png)

: **Figure 4**: Comparison of execution time and maximum resident set size (RSS) between `EIGENSTRAT` to `PACKEDANCESTRYMAP` (EGN-TO-PAM) and `PACKEDANCESTRYMAP`to `EIGENSTRAT` (PAM-TO-EGN) conversions for `dotgenolib`-based programs and `convertf` respectively.

While, in both cases, the software run time scales with the total number of elements (SNPs×individuals, N<sup>2</sup> in our benchmarks) in the file (**Figure 4**), there is a stark contrast in the case of peak memory usage. While `convertf` memory usage scales with the total number of elements in the file (N<sup>2</sup>), the `dotgenolib`-based methods scale with the number of _individuals_ (N). This can be see when we compare the relationship between the number of SNPs/individuals in a file and memory usage. While `convertf` scales polynomially, the `dotgenolib` based methods grow linearly (**Figure 4**). This difference is apparent when we compare peak memory usage when converting the largest `EIGENSTRAT` file (50,000 individuals and 50,000 SNPs) to the `PACKEDANCESTRYMAP` format. With `convertf` ~662,000 kB was required for this while only ~18,700 kB was required for the `dotgenolib`-based conversion (**Table 1**) - an over 35× decrease!

| Conversion | SNPs | `convertf` time (ms) | `convertf` max RSS (kB) | `dotgenolib` time (ms) | `dotgenolib` max RSS (kB) |
|---|---:|---:|---:|---:|---:|
| EGN → PAM | 100 | 0 | 2,628 | 0 | 1,876 |
| EGN → PAM | 500 | 20 | 3,008 | 0 | 1,904 |
| EGN → PAM | 1,000 | 50 | 3,696 | 10 | 2,076 |
| EGN → PAM | 2,000 | 210 | 5,428 | 50 | 2,312 |
| EGN → PAM | 3,000 | 410 | 7,792 | 50 | 2,572 |
| EGN → PAM | 5,000 | 1,390 | 13,516 | 130 | 3,336 |
| EGN → PAM | 8,000 | 5,260 | 26,040 | 360 | 4,820 |
| EGN → PAM | 10,000 | 7,430 | 36,780 | 540 | 5,120 |
| EGN → PAM | 15,000 | 20,030 | 72,168 | 1,150 | 7,092 |
| EGN → PAM | 30,000 | 125,680 | 251,676 | 4,530 | 12,940 |
| EGN → PAM | 40,000 | 276,830 | 433,020 | 8,300 | 15,820 |
| EGN → PAM | 50,000 | 668,300 | 662,256 | 12,230 | 18,704 |
| PAM → EGN | 100 | 0 | 2,616 | 170 | 1,672 |
| PAM → EGN | 500 | 30 | 2,964 | 10 | 2,008 |
| PAM → EGN | 1,000 | 40 | 3,824 | 30 | 1,920 |
| PAM → EGN | 2,000 | 150 | 5,404 | 40 | 2,316 |
| PAM → EGN | 3,000 | 270 | 7,736 | 80 | 2,688 |
| PAM → EGN | 5,000 | 810 | 13,900 | 210 | 3,416 |
| PAM → EGN | 8,000 | 1,990 | 26,412 | 560 | 4,388 |
| PAM → EGN | 10,000 | 3,370 | 36,800 | 860 | 5,188 |
| PAM → EGN | 15,000 | 7,810 | 71,680 | 1,800 | 7,180 |
| PAM → EGN | 30,000 | 38,810 | 250,884 | 7,080 | 12,888 |
| PAM → EGN | 40,000 | 90,050 | 431,396 | 12,760 | 15,628 |
| PAM → EGN | 50,000 | 173,970 | 660,428 | 19,700 | 18,572 |

: **Table 1**: Run time and memory usage benchmarks for `EIGENSTRAT` (EGN) and `PACKEDANCESTRYMAP` (PAM) conversions

This difference can be explained by the different approaches each piece of software takes. While `convertf` loads the entire genotype matrix into memory, `dotgenolib` allows one to easily use the simpler, more resourceful route of reading the file as a stream. Not only does this approach use less memory than `convertf`, it runs faster too!

### Merging inputs
While the AADR is a very useful dataset, it is (like all datasets) limited in scope. There are numerous populations, both ancient and modern, which are outside of its scope and, because of this, researchers who are interested in comparing their genomic datasets with those in the AADR have a non-trivial issue to overcome - merging datasets. As alluded to previously, the `PACKEDANCESTRYMAP` format is the _de facto_ standard for a wide variety of workflows in the world of population genetics, so datasets should ideally be able to be merged into this format. Thankfully, `AdmixTools` comes with a command for just this task! The `mergeit` command can merge two `PACKEDANCESTRYMAP` files into a single file with samples from both datasets. This software, however, comes with the same pitfalls as `convertf`. Much like `convertf`, it is VERY memory hungry, requiring that data from **both** files be loaded into memory. For large `PACKEDANCESTRYMAP` files like the AADR this basically guarantees that merging will be impossible on a standard desktop or laptop. In fact, this was one of **the** motivating factors behind me making this package.

In response to these issues I created a program, [`pam_merge`](https://github.com/floutt/admix_toolbox/blob/main/src/pam_merge.c), for merging multiple `PACKEDANCESTRYMAP` files. In benchmarks of speed and memory use it outperformed `mergeit`, running up to 4.9× faster and with an order of magnitude less peak memory usage (**Table 2**, **Figure 5**). Benchmarks were performed by merging v54 of the AADR (the Human Origins SNP set, `aadr_v54.1.p1_HO_public*`) with simulated `PACKEDANCESTRYMAP` files with 12 samples and variable numbers of SNPs and comparing run time and memory usage between `mergeit` and `pam_merge` with the `--intersect-snps` option.

![]({static}/images/dotgeno/mergeit_cmp.png)

: **Figure 5**: Comparison of execution time and peak memory usage between the `dotgenolib`-based `pam_merge` and `mergeit`.

Across the board `pam_merge` showed superior performance to `mergeit` (**Figure 5**), this difference is especially apparent in memory usage. While `pam_merge` was 2-4× faster than `mergeit`, it used 15-52× less memory (**Table 2**). These metrics are a night and day difference and allow people with more constrained hardware to work with datasets like the AADR.

| SNPs | `dotgenolib` time (ms) | `dotgenolib` max RSS (kB) | `mergeit` time (ms) | `mergeit` max RSS (kB) |
|---:|---:|---:|---:|---:|
| 500 | 950 | 124,812 | 4,670 | 6,477,296 |
| 1,000 | 960 | 125,500 | 4,620 | 6,477,060 |
| 2,000 | 1,110 | 125,928 | 5,010 | 6,476,952 |
| 3,000 | 1,210 | 126,400 | 5,270 | 6,477,160 |
| 5,000 | 1,390 | 127,800 | 5,590 | 6,476,908 |
| 8,000 | 1,700 | 129,868 | 6,570 | 6,477,284 |
| 10,000 | 1,960 | 130,940 | 7,350 | 6,477,108 |
| 15,000 | 2,390 | 134,816 | 7,760 | 6,477,052 |
| 30,000 | 4,010 | 145,116 | 10,440 | 6,477,124 |
| 40,000 | 5,050 | 150,988 | 13,310 | 6,477,136 |
| 50,000 | 6,170 | 156,740 | 16,130 | 6,476,984 |
| 100,000 | 11,290 | 188,704 | 29,170 | 6,477,164 |
| 500,000 | 53,960 | 463,156 | 113,500 | 6,752,040 |

: **Table 2**: Run time and memory usage benchmarks for merging v54 of the AADR with simulated `PACKEDANCESTRYMAP` files of varying sizes using the `dotgenolib`-based `pam_merge` and `mergeit`. 

### PLINK-style input filtering
Another big issue with both the `PACKEDANCESTRYMAP` and `EIGENSTRAT` formats is that there are few tools for performing [PLINK-style](https://www.cog-genomics.org/plink/2.0/filter) input filtering with them. In response to this, I developed [`admix_io`](https://github.com/floutt/admix_io) a piece of software which uses `dotgenolib` to do PLINK-style input filtering on this file type. To my surprise, parallel to my work on `dotgenolib`, [Dr. Dilek Koptekin](https://www.koptekin.org/) created [`pygenstrat`](https://github.com/dkoptekin/pygenstrat) for this very task! This software is written in Python and uses NumPy-cacheing to boost performance. One issue with this package is that it (as of the writing of this article) [does not support PACKEDANCESTRYMAP hashing](https://github.com/dkoptekin/pygenstrat/issues/3) which can be problematic for downstream software that requires this metadata in the header.

![]({static}/images/dotgeno/filt_cmp.png)

: **Figure 6**: Comparison of execution time and peak memory usage between `pygenstrat` and the `dotgenolib`-based `admix_io` software with individual (IND) and SNP filtering tasks for files of varying sizes. Individual and SNP counts are varied for the IND and SNP benchmarks respectively.

I performed benchmarks of standard input filtering tasks with both pieces of software. For this I generated simulated PACKEDANCESTRYMAP datasets with varying numbers of SNPs and individuals and ran tests for SNP and individual filtering respectively. For SNP filtering I filtered SNPs to only those in chromosome 1 with a minor allele frequency greater than or equal to 0.1 while for individual filtering samples from a single population were selected. For SNP filtering this was done on `PACKEDANCESTRYMAP` files with increasing numbers of SNPs (with each file having 5,000 individuals). Similarly, for individual filtering this was fone on files with increasing numbner of individuals (with each file having 5,000 SNPs). Run time and memory usage statistics were then compared between the `dotgenolib`-based `admix_io` and `pygenstrat`.

| Filter | Input size | `admix_io` time (ms) | `admix_io` max RSS (kB) | `pygenstrat` time (ms) | `pygenstrat` max RSS (kB) |
|---|---:|---:|---:|---:|---:|
| SNP | 5,000 | 40 | 4,528 | 1,850 | 148,340 |
| SNP | 10,000 | 50 | 5,604 | 1,360 | 152,036 |
| SNP | 50,000 | 280 | 14,072 | 3,200 | 162,408 |
| SNP | 70,000 | 150 | 17,712 | 3,790 | 161,068 |
| IND | 5,000 | 80 | 4,856 | 1,250 | 147,796 |
| IND | 10,000 | 130 | 5,344 | 1,460 | 178,196 |
| IND | 50,000 | 490 | 11,612 | 2,860 | 390,940 |
| IND | 70,000 | 630 | 15,848 | 3,700 | 500,760 |

: **Table 3**: Run time and memory usage benchmarks for an individual (IND) and SNP filtering task with the `dotgenolib`-based `admix_io` and `pygenstrat`

In these benchmarks `admix_io` universally showed superior performance to `pygenstrat` (**Figure 6**), with `admix_io` runs on the larger files showing superior performance to `pygenstrat` runs on the smallest files in regards to both speed and memory usage (**Table 3**). 

### Converting TGENO file to PACKEDANCESTRYMAP format
The AADR has recently made the move from distributing the AADR in the `PACKEDANCESTRYMAP` format to doing so in the `TGENO` one. This rapid change has added another hurdle to using the dataset, as more programs use the `PACKEDANCESTRYMAP` than the transposed `TGENO` one. As such, conversion to the latter format is a often a requirement. This change has also [exposed I/O bugs in existing programs](https://github.com/uqrmaie1/admixtools/issues/172) relating to `TGENO` files. Because of this, it is important to have a way to convert between `TGENO` and `PACKEDANCESTRYMAP` formats. While the `convertf` program allows one to perform this operation, it is very memory intensive (as noted earlier in this article). In response to this, `AdmixTools` now comes with a [specialized program](https://github.com/DReichLab/AdmixTools/blob/master/src/transpose.c) for converting `TGENO` files to the `PACKANCESTRYMAP` format, allowing users to either run the conversion by loading all of the data into memory or using the `--low-mem` flag to do this in a memory-thrifty manner. Both options come with a trade-offs, the fast method eats up memory (a precious resource for the hobbyist) while the `--low-mem` method is incredibly slow. As is often the case in life, the best tradeoff is usually somewhere in the middle, an option unavailable to the user of this program. In response to this, I used `dotgenolib` to develop a [program](https://github.com/floutt/admix_toolbox/blob/main/src/transpose_tgeno.c) to convert `TGENO` files to the `PACKEDANCESTRYMAP` format. This program allows the user to determine the number of SNP vectors to store in memory, thus allowing the user to determine what their own personal happy medium is.

### Converting VCF files to the PACKEDANCESTRYMAP format (and vice versa)
One of the most widely used genomic file formats is the Variant Call Format (VCF). This format is the standard for storing data obtained from whole-genome sequencing datasets. For example, it is used by projects such as the [1000 Genomes Project](https://www.internationalgenome.org/home/) and the [Genome Aggregation Database](https://gnomad.broadinstitute.org/) to store the genomic data of thousands of individuals. As such, many researchers and hobbyists may find it useful to convert VCF files into the `PACKEDANCESTRYMAP` format to ultimately merge this data with AADR for a variety of downstream analyses. Similarly, VCF files are required input for many pieces of software, making them useful for researchers who would like to perform particular analyses on samples from the AADR. For both of these reasons, having a method for the quick and easy conversion between the two formats is incredibly useful. In response to this I created a program, [`vcf_pam_conv`](https://github.com/floutt/admix_toolbox/blob/main/src/vcf_pam_conv.c), for converting `PACKEDANCESTRYMAP` files to VCF and VCF files to the `PACKEDANCESTRYMAP` format using `dotgenolib`.

## Concluding remarks
`dotgenolib` is a useful library for reading, writing, and manipulating the `PACKEDANCESTRYMAP`, `TGENO`, and `EIGENSTRAT` files. It has shown its utility in a variety of arenas from [allowing PLINK-style input filtering of `PACKEDANCESTRYMAP` and `EIGENSTRAT` files](https://github.com/floutt/admix_io), to providing [memory efficient alternatives](https://github.com/floutt/admix_toolbox) to `mergeit` and `convertf`. In all of these examples, the `dotgenolib`-based implementations showed superior performance in terms of memory usage and speed to existing methods. I hope that researchers and hobbyists alike find `dotgenolib`, and the tools written in it, to be useful!

## References
1. Eren, K., Taktakoğlu, N. & Pirim, I. DNA Sequencing Methods: From Past to Present. Eurasian J Med 54, S47–S56 (2022).
2. Lenoir, T. & Giannella, E. The emergence and diffusion of DNA microarray technology. J Biomed Discov Collab 1, 11 (2006).
3. Hagelberg, E., Hofreiter, M. & Keyser, C. Ancient DNA: the first three decades. Philos Trans R Soc Lond B Biol Sci 370, 20130371 (2015).
4. Yen, D. Sweet-potato variation and its relation to human migration in the Pacific. in 24 (Bishop Museum Press, 1963).
5. Hather, J. & Kirch, P. V. Prehistoric sweet potato (Ipomoea batatas) from Mangaia Island, Central Polynesia. Antiquity 65, 887–893 (1991).
6. Dunis, S. Of Kumara and Canoes: Maori and Hawaiian Mythologies and American Contacts. in The Sweet Potato in Oceania: A Reappraisal (eds Ballard, C., Brown, P., Bourke, R. M. & Harwood, T.) 89–97 (University of Sydney, Sydney, Australia, 2005).
7. Scaglion, R. Kumara in the Ecuadorian Gulf of Guayaquil? in The Sweet Potato in Oceania: A Reappraisal (eds Ballard, C., Brown, P., Bourke, R. M. & Harwood, T.) 35–41 (University of Sydney Press, Pittsburgh and Sydney, 2005).
8. Ioannidis, A. G. et al. Native American gene flow into Polynesia predating Easter Island settlement. Nature 583, 572–577 (2020).
9. Moreno-Mayar, J. V. et al. Ancient Rapanui genomes reveal resilience and pre-European contact with the Americas. Nature 633, 389–397 (2024).
10. Mallick, S. et al. The Allen Ancient DNA Resource (AADR) a curated compendium of ancient human genomes. Sci Data 11, 182 (2024).

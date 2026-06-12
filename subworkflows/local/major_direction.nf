include {map_to_build} from '../../modules/local/map_to_build'
include {ten_percent_counts} from '../../modules/local/ten_percent_counts'
include {ten_percent_counts_sum} from '../../modules/local/ten_percent_counts_sum'
include {generate_strand_counts} from '../../modules/local/generate_strand_counts'
include {summarise_strand_counts} from '../../modules/local/summarise_strand_counts'

workflow major_direction{
    take:
    chr
    files

    main:
    // Fan out: each (sample, chromosome) pair runs its own map_to_build process.
    // This replaces the single all-chromosome job, reducing per-process memory ~7x
    // and allowing all chromosomes to map in parallel.
    chrom_ch = chr.flatten().map { it.toString().replaceAll("chr", "") }
    map_to_build(files.combine(chrom_ch))

    // Reshape output: [GCST, chrom, merged, unmapped, yaml]
    //   → map_chr_ch: [chrN, GCST, merged, yaml]  (same key layout as before)
    map_to_build.out.mapped
        .map { gcst, chrom, merged, unmapped_file, yaml ->
            tuple("chr" + chrom, gcst, merged, yaml)
        }
        .set { map_chr_ch }

    // Merge per-chromosome unmapped files into one file per GCST for the log.
    // collectFile(keepHeader:true) concatenates content, keeping the header once.
    unmapped = map_to_build.out.mapped
        .map { gcst, chrom, merged, unmapped_file, yaml ->
            tuple(gcst, unmapped_file)
        }
        .collectFile(keepHeader: true) { gcst, uf ->
            ["${gcst}.unmapped", uf.text]
        }
        .map { f -> tuple(f.getBaseName().replace('.unmapped', ''), f) }

    Channel.fromPath("${params.ref}/homo_sapiens-chr*.vcf.gz")
           .map { prepare_reference(it) }
           .set { ref_chr_ch }

    count_ch = map_chr_ch.combine(ref_chr_ch, by: 0)

    ten_percent_counts(count_ch)

    int nchr = params.chrom.size()
    ten_to_sum = ten_percent_counts.out
                      .ten_sc
                      .groupTuple(by: 0)
                      .branch { pass: it[1].size() == nchr }
                      .map { it[0] }

    ten_percent_counts_sum(ten_to_sum)

    ten_percent_counts_sum.out.ten_sum.branch { rerun:   it.contains("rerun")
                                                contiune: it.contains("contiune") }
                                      .set { branch }

    branch.rerun.map { tuple(it[3], it[0]) }.set { all_sc_ch }
    count_ch.combine(all_sc_ch, by: 1).set { rerun_branch }
    generate_strand_counts(rerun_branch)

    all_to_sum = generate_strand_counts.out.all_sc.collect().map { tuple(it[0], it[1]) }.unique()
    summarise_strand_counts(all_to_sum)

    all_files = summarise_strand_counts.out.all_sum.mix(branch.contiune)

    rearrnaged_count_ch = count_ch.map { tuple(it[1], it[0], it[2], it[3], it[4]) }
    all_input = all_files.combine(rearrnaged_count_ch, by: 0)
    hm_input = all_input.map { it[0, 2..7] }
    direction_sum = all_input.map { it[0..1] }.unique()

    emit:
    hm_input     = hm_input
    direction_sum = direction_sum
    unmapped     = unmapped
}

// groovy helpers
def prepare_reference(Path input) {
    return [input.getName().split('-')[1].split('\\.')[0], input]
}

def get_chr(Path input) {
    return ("chr" + input.getName().split('\\.')[0])
}

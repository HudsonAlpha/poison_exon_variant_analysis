workflow bcftools_filter_bedtools_intersect {
  meta {
      author: "Stephanie Felker"
      email: "sfelker@hudsonalpha.org"
      description: "Filter a list of vcfs for variants within a input bed, removes headers, and replaces them with a dummy header for downstream analysis."
  }

  input {
      Array[File] VCF_LIST
      String COHORT
      File REGION_BED
      Int CORES = 1
      Int DISK = 100
      Int MEM = 10
      Int PREEMPTIBLE = 0
  }

  call filter_intersect {
      input:
      in_vcf_list=VCF_LIST,
      cohort=COHORT,
      region_bed=REGION_BED,
      in_cores=CORES,
      in_disk=DISK,
      in_mem=MEM,
      in_preemptible=PREEMPTIBLE
  }
}

task filter_intersect {
    input {
      Array[File] in_vcf_list
      File region_bed
      String cohort
      Int in_cores
      Int in_disk
      Int in_mem
      Int in_preemptible
    }
    command <<<
    set -eux -o pipefail
    # prepare scripts
        rm -f vcf_list.txt
        mkdir -p temp
        for x in ~{sep=' ' in_vcf_list}
        do
            # QC filtering
            bcftools filter -i 'FILTER="PASS"' $x | bcftools filter -i 'FORMAT/DP>10' > $(basename "${x}").initfilter.vcf
         #   # Intersect with regions of interest
            bedtools intersect -a $(basename "${x}").initfilter.vcf -b ~{region_bed} -header -wa | uniq | bgzip > $(basename "${x}").regionfilter.vcf.gz
            # index files
            tabix -p vcf $(basename "${x}").regionfilter.vcf.gz
            # remove annotations and change headers
            bcftools annotate -x INFO,FORMAT $(basename "${x}").regionfilter.vcf.gz | sed "/^##/d"| sed "1s/^/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n/" | bgzip > $(basename "${x}").final.vcf.gz
            # clean up and index
            tabix -p vcf $(basename "${x}").final.vcf.gz
           # rm -f $(basename "${x}").initfilter.vcf $(basename "${x}").regionfilter.vcf.gz $(basename "${x}").regionfilter.vcf.gz.tbi
        done;

    >>>
    output {
        Array[File] out_vcf= glob("*.final.vcf.gz")
        Array[File] out_vcf_idx= glob("*.final.vcf.gz.tbi")


    }
    runtime {
        docker: "safelker/richard_cser_anvil:21_Dec_21_update"
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        preemptible: in_preemptible
    }
}

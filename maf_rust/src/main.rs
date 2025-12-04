use std::env;
use std::process;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use rust_htslib::{bam, bam::Read};

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn main() {
  let filter_1kg = true;

  let t2t_chroms:[&str;24] = ["NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1", "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1", "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1", "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1", "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1"];
  let chroms:[&str;24] = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"];

  let args: Vec<String> = env::args().collect();
  if args.len() < 5 {
    eprintln!("usage: nano_maf <aligned_sorted.BAM> <enriched_regions.BED> <repetitive_regions.BED> <variable_sites.TSV> <min_depth>");
    process::exit(1);
  }

  let mut regions = Vec::new();
  for _c in 0..24 {
      regions.push(Vec::new());
  }
  if let Ok(lines) = read_lines(&args[2]) {
    // Consumes the iterator, returns an (Optional) String
    for line in lines.flatten() {
      let parts:Vec<&str> = line.trim().split("\t").collect();
      let chr = parts[0];
      let index = {
        if t2t_chroms.contains(&chr) {
          t2t_chroms.iter().position(|&r| r == chr).unwrap()
        } else if chroms.contains(&chr) {
          chroms.iter().position(|&r| r == chr).unwrap()
        } else {
          eprintln!("WARNING: Unknown chromosome in repeats file (skipping it): {}", chr);
          continue;
        }
      };
      let st:u32 = parts[1].parse().unwrap();
      let en:u32 = parts[2].parse().unwrap();
      if (parts[3].len() >= 3 && parts[3][..3] == *"IGH") || (parts[3].len() >=4 && parts[3][..4] == *"DUX4") { // hard ignore these complex regions that do not behave for variant calling
        continue;
      }
      //println!("{}/{}: {}", chr, index, pos);
      regions[index].push((st, en));
    }
  } else {
    eprintln!("ERROR: failed to read BED file");
    process::exit(1);
  }

  // read repeat BED
  let mut repeats = Vec::new();
  for _c in 0..24 {
      repeats.push(Vec::new());
  }
  if let Ok(lines) = read_lines(&args[3]) {
    for line in lines.flatten() {
      let parts:Vec<&str> = line.trim().split("\t").collect();
      let chr = parts[0];
      let index = {
        if t2t_chroms.contains(&chr) {
          t2t_chroms.iter().position(|&r| r == chr).unwrap()
        } else if chroms.contains(&chr) {
          chroms.iter().position(|&r| r == chr).unwrap()
        } else {
          eprintln!("WARNING: Unknown chromosome in repeats file (skipping it): {}", chr);
          continue;
        }
      };
      let st:u32 = parts[1].parse().unwrap();
      let en:u32 = parts[2].parse().unwrap();
      repeats[index].push((st, en));
    }
  } else {
    eprintln!("ERROR: failed to read repeat BED file");
    process::exit(1);
  }
  for i in 0..24 {
    repeats[i].sort();
  }

  // read 1Kgenomes sites (not BED, but close)
  let mut sites = Vec::new();
  for _c in 0..24 {
      sites.push(Vec::new());
  }
  if let Ok(lines) = read_lines(&args[4]) {
    for line in lines.flatten() {
      let parts:Vec<&str> = line.trim().split("\t").collect();
      let chr = parts[0];
      let index = {
        if t2t_chroms.contains(&chr) {
          t2t_chroms.iter().position(|&r| r == chr).unwrap()
        } else if chroms.contains(&chr) {
          chroms.iter().position(|&r| r == chr).unwrap()
        } else {
          eprintln!("WARNING: Unknown chromosome in sites file (skipping it): {}", chr);
          continue;
        }
      };
      if parts.len() < 4 {
        eprintln!("ERROR parsing sites file: {:#?}", parts);
      }
      let pos:u32 = parts[1].parse().unwrap();
      let al0:u32 = match parts[2] {
        "A" => 0,
        "C" => 1,
        "G" => 2,
        "T" => 3,
        _ => 4
      };
      let al1:u32 = match parts[3] {
        "A" => 0,
        "C" => 1,
        "G" => 2,
        "T" => 3,
        _ => 4
      };
      sites[index].push((pos, al0, al1));
    }
  } else {
    eprintln!("ERROR: failed to read repeat BED file");
    process::exit(1);
  }

  let min_depth = match args[5].parse::<i32>() {
      Ok(num) => num,
      Err(e) => panic!("min_depth cannot be interpreted as an integer: {e}"),
  };

  let mut bam = bam::IndexedReader::from_path(&args[1]).unwrap();

  eprintln!("{}", &args[1]);
  for chrom in chroms.iter() {
    let chr_idx = chroms.iter().position(|r| r == chrom).unwrap();
    let ref_names = bam.header().target_names();
    let chrom_name = if ref_names.into_iter().find(|&a| a==chrom.as_bytes()) != None {
        chrom
    } else {
        t2t_chroms[chr_idx]
    };
    eprintln!("  {chrom_name} ({chrom}) ...");

    let mut g:usize = 0; // index into regions
    let mut s:usize = 0; // index into sites
    let mut r:usize = 0; // index into repeats

    // iterate over BED regions in chrom
    while g < regions[chr_idx].len() {
      eprintln!("    region {g} ({} - {}) ...", regions[chr_idx][g].0, regions[chr_idx][g].1);

      let _ = bam.fetch((chrom_name, regions[chr_idx][g].0, regions[chr_idx][g].1)).expect("Couldn't fetch from BAM file");
      // pileup over all covered sites
      for p in bam.pileup() {
        let pileup = p.unwrap();
        let pos = pileup.pos()+1; // change from 0-based to 1-based index, to match coninical genome annotation
        if pos < regions[chr_idx][g].0 {
          continue;
        }
        if pos > regions[chr_idx][g].1 {
          break;
        }
        //eprintln!("{}:{} depth {}", pileup.tid(), pos, pileup.depth());

        // move repeats index (r) up to current or next repeat entry
        while r < repeats[chr_idx].len() && repeats[chr_idx][r].1 < pos {
            r += 1;
        }
        // skip site if in a repeat element
        if pos >= repeats[chr_idx][r].0 && pos <= repeats[chr_idx][r].1 {
            //eprintln!("Skip {chrom} @ {pos}");
            continue;
        }

        if filter_1kg {
            // move sites index (s) up to current or next site entry
            while s < sites[chr_idx].len() && sites[chr_idx][s].0 < pos {
                s += 1;
            }
            // skip site if NOT a 1KG variable site
            if s >= sites[chr_idx].len() || pos != sites[chr_idx][s].0 {
                continue;
            }
        }
    
        let mut allele_counts: [i32; 5] = [0; 5]; // A, C, G, T, N/? counts
        for alignment in pileup.alignments() {
          if alignment.record().is_secondary() { // secondary alignments do not include the read sequence
              continue;
          }
          if !alignment.is_del() && !alignment.is_refskip() {
            //eprintln!("Seq: {}", alignment.record().seq().as_bytes().len());
            let nt:char = alignment.record().seq()[alignment.qpos().unwrap()].try_into().unwrap();
            let al = match nt {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => 4
            };
            allele_counts[al] += 1;
          }
          /*
          match alignment.indel() {
            bam::pileup::Indel::Ins(len) => eprintln!("Insertion of length {} between this and next position.", len),
            bam::pileup::Indel::Del(len) => eprintln!("Deletion of length {} between this and next position.", len),
            bam::pileup::Indel::None => ()
          }
          */
        }
        let mut maj = 1; // highest count allele
        let mut min = 0; // SECOND highest count allele
        if allele_counts[min] > allele_counts[maj] {
            min = 1;
            maj = 0;
        }
        for i in 2..4 {
            if allele_counts[i] > allele_counts[maj] {
                min = maj;
                maj = i
            } else if allele_counts[i] > allele_counts[min] {
                min = i;
            }
        }
        //eprintln!("min {}, maj {}", allele_counts[min], allele_counts[maj]);
        let depth = allele_counts[maj] + allele_counts[min];

        // Ref allele matches 1KG annotation, sometimes for errors the alt allele does not match,
        // these should be ignored
        //eprintln!("{chrom} {pos}: {} -> {}, maj: {}, min: {}", sites[chr_idx][s].1, sites[chr_idx][s].2, maj, min);
        if filter_1kg {
            if (maj == sites[chr_idx][s].1 as usize && min != sites[chr_idx][s].2 as usize) || (maj == sites[chr_idx][s].2 as usize && min != sites[chr_idx][s].1 as usize) {
                continue;
            }
        }

        // skip problematic section of chrom 5 (T2T)
        if chr_idx == 4 && pos > 720000 && pos < 760000 {
          continue;
        }

        if depth >= min_depth {
          println!("{chrom_name}\t{pos}\t{}\t{}", allele_counts[maj], allele_counts[min]);
        }
      }
      g += 1;
    }

  }
}

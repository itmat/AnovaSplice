usage= <<-eof
USAGE:
  ruby fix_gtf.rb gtf_file mapping_file fixed_gtf_file
  ____________________________________________________
  Fixes gene_id in gtf file
  IN: gtf_file mapping_file
  OUT: fixed_gtf_file
  mapping_file needs to look like this:
  #Format: tax_id GeneID Ensembl_gene_identifier RNA_nucleotide_accession.version Ensembl_rna_identifier protein_accession.version Ensembl_protein_identifier (tab is used as a separator, pound sign - start of a comment)
  7227  30970 FBgn0040373 NM_001297803.1  FBtr0340207 NP_001284732.1  FBpp0309182
  7227  30970 FBgn0040373 NM_130477.4 FBtr0070108 NP_569833.1 FBpp0070103
  7227  30970 FBgn0040373 NM_166834.2 FBtr0070107 NP_726658.1 FBpp0070102
  ____________________________________________________

eof

if ARGV.length != 3
  puts usage
  exit(2)
end

mapping = {}
File.open(ARGV[1]).each do |line|
  line.chomp!
  line = line.split(" ")
  mapping[line[3].split(".")[0]] = line[2]
end

#puts mapping["NM_166834"] #should be  FBgn0040373
#exit

out_file = File.open(ARGV[2],'w')

File.open(ARGV[0]).each do |line|
  line.chomp!
  #puts line
  line =~ /(NM_\d*)/
  id = $1
  if mapping[id]
    out_file.puts line.sub(/NM_\d*/,mapping[id])
  end
end

out_file.close()

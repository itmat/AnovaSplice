usage= <<-eof
USAGE:
  ruby fix_gtf.rb gtf_file mapping_file fixed_gtf_file
  ____________________________________________________
  Fixes gene_id in gtf file
  IN: gtf_file mapping_file
  OUT: fixed_gtf_file
  mapping_file needs to look like this:
   #name name2
   ENSMUST00000072177  ENSMUSG00000009772
   ENSMUST00000082125  ENSMUSG00000009772
   ENSMUST00000132064  ENSMUSG00000025909
   ENSMUST00000140295  ENSMUSG00000025909
   ENSMUST00000140302  ENSMUSG00000025909
   ENSMUST00000115484  ENSMUSG00000025909
   ENSMUST00000135046  ENSMUSG00000025909
   ENSMUST00000027042  ENSMUSG00000025909
   ENSMUST00000115488  ENSMUSG00000025909
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
  mapping[line[0]] = line[1]
end

#puts mapping["ENSMUST00000027042"] # Should be ENSMUSG00000025909

out_file = File.open(ARGV[2],'w')

File.open(ARGV[0]).each do |line|
  line.chomp!
  #puts line
  line =~ /(ENSMUST\d*)/
  #puts $1
  out_file.puts line.sub(/ENSMUST\d*/,mapping[$1])
end

out_file.close()

require "csv"

usage= <<-eof
USAGE:
  ruby mann_whitney_to_track.rb.csv mann_whitney.csv out.bed
  ____________________________________________________
  Makes UCSC track for significant exons
  IN: mann_whitney.csv
  OUT: out.bed
  # mann whitney
  mann_whitney.csv needs to look like this:
    "","id","Sample_9574_SS","Sample_9575_SS","Sample_9578_SS","Sample_9579_SS","Sample_9580_SS","Sample_9581_SS","Sample_9582_SS","Sample_9583_SS","P_Val","FDR"
    "1","exon:chr1:3204563-3207049",1,0,1,2,2,2,1,3,0.648941813187414,0.888219437653876
    "10","exon:chr1:4481009-4482749",24,24,20,10,11,16,13,19,0.46782507728497,0.761472006238808
    "11","exon:chr1:4483181-4483547",4,4,5,1,4,2,1,2,0.36880340134896,0.68167550231776
    "12","exon:chr1:4483181-4483571",4,4,6,1,4,2,1,2,0.36880340134896,0.68167550231776
  ____________________________________________________

eof

if ARGV.length != 2
  puts usage
  exit(2)
end

cutoff = 0.25

out_bed = File.open(ARGV[1], "w")
out_bed.puts "track name=\"Signinficant Exons\" description=\"Output from ttest_to_track.rb\""
CSV.foreach(ARGV[0], headers: true) do |row|
  next unless row["FDR"].to_f < cutoff
  chrom = row["id"].split(":")[1]
  start = row["id"].split(":")[2].split("-")[0].to_i
  stop = row["id"].split(":")[2].split("-")[1].to_i
  #chrom chromStart chromEnd name score(0-1000) strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
  # chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
  out_bed.puts "#{chrom}\t#{start}\t#{stop}"
end

out_bed.close
require "csv"

usage= <<-eof
USAGE:
  ruby sign_exons.rb anno.gtf mann_whitney.csv experiment.csv results.csv sign_exons.csv
  ____________________________________________________
  Adds number of exons to in.csv for ENSMUSG
  IN: file.gtf in.csv
  OUT: out.csv
  head -5 results.csv
    "","X","ensgene","pExon","pGroup","pInteraction","NumberOfExons","NumberOfSpliceforms","pExonFDR","pGroupFDR","pInteractionFDR"
    "1",1,"ENSMUSG00000009772",0.153753217474881,0.165372779024116,0.844479907905903,8,2,0.182107461323043,0.283061253685133,1
    "2",2,"ENSMUSG00000047528",0.0118258017270025,0.00989144699745502,0.187842966842855,15,1,0.0155413953091948,0.0295472184738957,0.752544141458701
    "3",3,"ENSMUSG00000026134",3.64302659836332e-10,0.273832557832243,0.90811608560453,14,1,7.72438527406982e-10,0.41094697897545,1
    "4",4,"ENSMUSG00000033569",NA,NA,NA,31,1,NA,NA,NA

  head -5 mann_whitney.csv
    "","id","Sample_9574_SS","Sample_9575_SS","Sample_9578_SS","Sample_9579_SS","Sample_9580_SS","Sample_9581_SS","Sample_9582_SS","Sample_9583_SS","P_Val","FDR"
    "1","exon:chr1:3204563-3207049",1,0,1,2,2,2,1,3,0.648941813187414,0.888219437653876
    "10","exon:chr1:4481009-4482749",24,24,20,10,11,16,13,19,0.46782507728497,0.761472006238808
    "11","exon:chr1:4483181-4483547",4,4,5,1,4,2,1,2,0.36880340134896,0.68167550231776
    "12","exon:chr1:4483181-4483571",4,4,6,1,4,2,1,2,0.36880340134896,0.68167550231776
  ____________________________________________________

eof

cutoff = 0.25

if ARGV.length != 5
  puts usage
  exit(2)
end

## Read anno.gtf
mapping = {}
File.open(ARGV[0]).each do |line|
  line.chomp!
  line = line.split("\t")
  next unless line[2] == "exon"
  line[-1] =~ /(ENSMUSG\d*)/
  name = $1.delete("\"")
  mapping["exon:#{line[0]}:#{line[3]}-#{line[4]}"] ||= []
  mapping["exon:#{line[0]}:#{line[3]}-#{line[4]}"] << name unless mapping["exon:#{line[0]}:#{line[3]}-#{line[4]}"].include?(name)
end

#puts mapping["exon:chr1:4481009-4482749"] # should be ["ENSMUSG00000025902"]
#exit


# Read experiment.csv
experiment = {}
CSV.foreach(ARGV[2], headers: true) do |row|
  experiment[row["Group"]] ||= []
  experiment[row["Group"]] << row["SampleName"]
end


# Read ttest.csv
counter = 0
current_gene = nil
ens_gene = nil
keeper = []
upreg = false
downreg = false
exons_up = {}
exons_down = {}
CSV.foreach(ARGV[1], headers: true) do |row|
  next unless row["FDR"].to_f < cutoff
  #puts row["FDR"]
  next unless mapping[row["id"]]
  ens_gene = mapping[row["id"]][0]
  current_gene ||= ens_gene
  if ens_gene != current_gene
    if upreg && downreg
      keeper << current_gene
      counter += 1
    else
      exons_up.delete(current_gene)
      exons_down.delete(current_gene)
    end
    current_gene = ens_gene
    upreg = false
    downreg = false
  end

  sum = Array.new(experiment.keys.length,0)
  i = 0
  experiment.keys.each do |g|
    experiment[g].each do |sample_name|
      sum[i] += row[sample_name].to_i
    end
    i += 1
  end
  #puts sum.join(":")
  if sum[0] > sum[1]
    exons_up[current_gene] ||= []
    exons_up[current_gene] << row["id"]
    upreg = true
  else #sum[0] < sum[1]
    exons_down[current_gene] ||= []
    exons_down[current_gene] << row["id"]
    downreg = true
  end
end

if upreg && downreg
  keeper << ens_gene
  counter += 1
end

puts counter
puts keeper.join("\t")
puts keeper.include?("ENSMUSG00000025231")
puts exons_up
puts exons_down

#exit

CSV.open(ARGV[-1],'w') do |csv|
  csv << ["nana","nana",  "ensgene","pExon","pGroup","pInteraction","NumberOfExons","NumberOfSpliceforms","pExonFDR","pGroupFDR","pInteractionFDR"]
  CSV.foreach(ARGV[-2], headers: true) do |row|
    #puts "#{row["ensgene"]}"
    if keeper.include?(row["ensgene"])
      #puts "YESSSSSSS"
      csv.puts row
    end
  end
end

#out_file.close()

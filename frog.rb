#!/usr/bin/ruby

#There was a green grasshopper|3 p
#As a cucumber so green.

#Imagine youself|2 p
#There was a green grasshopper
#Imagine youself|2 p
#As a cucumber so green.

#It ate the only green grass|3 p
#And made good friends with flies.

#Imagine youself|2 p
#It ate the only green grass
#Imagine youself|2 p
#And made good friends with flies.

#But one day came a big frog |3 p
#And ate the hopper ! Ah!!!

#Imagine youself|2 p
#But one day came a big frog
#Imagine youself|2 p
#And ate the hopper ! Ah!!!

require "fileutils"

def get_opt(option, with_value = false)
  res = ARGV.index(option)
  return nil unless res
  return ARGV[res] unless with_value

  return nil if ARGV.length <= res + 1

  ARGV[res + 1]
end

modes = ["ECB", "CBC", "CFB", "OFB"]

if get_opt "--help"
  puts "-r N                   repeat each test N times"
  puts "-m \"mode1 mode2 ...\"   modes for testing. #{modes.join(', ')} by default"
  puts "-i file                specify input file. \"input.txt\" by default"
  puts "-q                     report only warnings and the result"
  exit 1
end

repeat = get_opt "-r", true
repeat = repeat ? repeat.to_i : 1

input = get_opt "-i", true
input = input ? input : "input.txt"

part_modes = get_opt "-m", true
modes = part_modes ? part_modes.split : modes

quiet = get_opt "-q"

dirname = "TEST_" + Time.now.asctime.gsub(/\W/, '-')
Dir.mkdir dirname

time = {}

Dir.chdir(dirname) do

  inp_file = "input.txt"
  FileUtils.copy "../grasshopper", "."
  FileUtils.copy "../#{input}", inp_file

  repeat.times do |iteration|

    puts "\nIteration #{iteration}" unless quiet

    modes.each do |opt|
      unless time[opt]
        time[opt] = []
      end

      exe_res = `./grasshopper -#{opt}`

      enc_file = "e#{opt}#{iteration}"
      dec_file = "d#{opt}#{iteration}"
      FileUtils.move "encrypted.txt", enc_file
      FileUtils.move "decrypted.txt", dec_file
      if not FileUtils.cmp inp_file, dec_file
        puts "\t#{opt}: #{inp_file} and #{dec_file} files are different"
      end

      exe_res = exe_res.lines.select{ |line| line =~ /Speed/ }
      exe_res.map!{ |line| line = line[/(\d+([.]\d*)?|int)/].to_f }
      exe_res.map!{ |time| time == 0 ? Float::INFINITY : time }

      enc_time = exe_res[0]
      dec_time = exe_res[1]

      puts "\t#{opt} encrypt speed: #{enc_time} MB/s" unless quiet
      puts "\t#{opt} decrypt speed: #{dec_time} MB/s" unless quiet

      if enc_time != Float::INFINITY && dec_time != Float::INFINITY
        time[opt] << {:enc_time => enc_time, :dec_time => dec_time}
      else
        puts "\tDon't add infinite #{opt} time result"
      end
      puts "" unless quiet
    end
  end

  puts "Average"
  m = "Mode"
  e = "Encrypt speed MB/s"
  d = "Decrypt speed MB/s"
  puts "| #{m} | #{e} | #{d} "
  time.each do |key, val|
    next if val.empty?
    enc_time = val.inject(0){ |sum, h| sum += h[:enc_time] }
    dec_time = val.inject(0){ |sum, h| sum += h[:dec_time] }

    enc_time = enc_time / val.size
    dec_time = dec_time / val.size

    puts "| #{key.ljust(m.size)} | #{enc_time.to_s.ljust(e.size)} | #{dec_time.to_s.ljust(d.size)}"
  end
end




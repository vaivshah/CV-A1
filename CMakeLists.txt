cmake_minimum_required(VERSION 3.9)
project(ramchand_a1)

set(CMAKE_CXX_STANDARD 11)

add_executable(ramchand_a1
        part1/denoised.png
        part1/DTwoDimArray.h
        part1/fft.h
        part1/Makefile
        part1/noise1.png
        part1/SImage.h
        part1/SImageIO.h
        part1/watermark
        part1/watermark.cpp
        part1/watermarked.png
        part1/wmspectro.png
        part2/detect
        part2/detect.cpp
        part2/detected.png
        part2/detected.txt
        part2/DTwoDimArray.h
        part2/ic_1.png
        part2/ic_2.png
        part2/ic_3.png
        part2/ic_4.png
        part2/ic_5.png
        part2/ic_6.png
        part2/ic_7.png
        part2/ic_8.png
        part2/ic_9.png
        part2/Makefile
        part2/SImage.h
        part2/SImageIO.h
        README.md)

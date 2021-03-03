// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <fstream>

#include "catch.hpp"
#include "helpers.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/files/GzFile.hpp"
#include "chemfiles/Error.hpp"
using namespace chemfiles;

static void check_file(TextFile& file) {
    CHECK(file.readline() == "297");
    CHECK(file.readline() == " generated by VMD");
    CHECK(file.readline() == "  O          0.417219        8.303366       11.737172");

    file.rewind();
    CHECK(file.readline() == "297");
    CHECK(file.readline() == " generated by VMD");

    // Count lines
    file.rewind();
    size_t lines = 0;
    while (!file.eof()) {
        file.readline();
        lines++;
    }

    CHECK(lines == 29901);
    CHECK(file.tellpos() == 1606000);
    CHECK(file.eof());

    file.seekpos(23804);
    CHECK(file.readline() == "  H          8.479585        0.521128       11.514298");
}

TEST_CASE("Read a text file") {
    SECTION("Read at different compressions levels") {
        auto file_6 = TextFile("data/xyz/water.6.xyz.gz", File::READ, File::GZIP);
        check_file(file_6);

        auto file_9 = TextFile("data/xyz/water.9.xyz.gz", File::READ, File::GZIP);
        check_file(file_9);
    }

    SECTION("Constructor errors") {
        CHECK_THROWS_WITH(
            GzFile("not existing", File::READ),
            "could not open the file at 'not existing'"
        );
    }

    SECTION("Lines offsets") {
        // Compare offset with uncompressed file
        auto file = TextFile("data/xyz/water.xyz", File::READ, File::DEFAULT);
        auto positions = std::vector<uint64_t>();
        while (!file.eof()) {
            positions.push_back(file.tellpos());
            file.readline();
        }

        auto gz_file = TextFile("data/xyz/water.6.xyz.gz", File::READ, File::GZIP);
        for (size_t i=0; i<positions.size(); i++) {
            CHECK(positions[i] == gz_file.tellpos());
            gz_file.readline();
        }
        CHECK(gz_file.eof());
    }
}

TEST_CASE("Write a gz file") {
    auto filename = NamedTempPath(".gz");

    {
        TextFile file(filename, File::WRITE, File::GZIP);
        file.print("Test\n");
        file.print("{}\n", 5467);
        CHECK(file.tellpos() == 10);
    }

    std::ifstream checking(filename, std::ios::binary);
    REQUIRE(checking.is_open());
    checking.seekg(0, std::ios::end);
    auto size = static_cast<size_t>(checking.tellg());
    checking.seekg(0, std::ios::beg);

    auto content = std::vector<uint8_t>(size);
    checking.read(reinterpret_cast<char*>(content.data()), static_cast<std::streamsize>(size));

    // Byte 9 identify the OS in gzip files.
    // Override it so we can check for the full file output
    content[9] = 0xff;

    auto expected = std::vector<uint8_t> {
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x0b, 0x49,
        0x2d, 0x2e, 0xe1, 0x32, 0x35, 0x31, 0x33, 0xe7, 0x02, 0x00, 0x8a, 0x43,
        0x5e, 0x98, 0x0a, 0x00, 0x00, 0x00,
    };
    CHECK(content == expected);

    // Decompress and compare
    TextFile file(filename, File::READ, File::GZIP);
    CHECK(file.readline() == "Test");
    CHECK(file.readline() == "5467");
}

TEST_CASE("Append to a gz file") {
    auto filename = NamedTempPath(".gz");

    {
        TextFile file(filename, File::APPEND, File::GZIP);
        file.print("Append 1\n");
        file.print("{}\n", 7645);
    }

    {
        TextFile file(filename, File::APPEND, File::GZIP);
        file.print("Append 2\n");
        file.print("{}\n", 6754);
    }

    std::ifstream checking(filename, std::ios::binary);
    REQUIRE(checking.is_open());
    checking.seekg(0, std::ios::end);
    auto size = static_cast<size_t>(checking.tellg());
    checking.seekg(0, std::ios::beg);

    auto content = std::vector<uint8_t>(size);
    checking.read(reinterpret_cast<char*>(content.data()), static_cast<std::streamsize>(size));

    // Byte 9 identify the OS in gzip files.
    // Override it so we can check for the full file output
    content[9] = 0xff;
    // A new header is written when appending
    content[43] = 0xff;

    auto expected = std::vector<uint8_t> {
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
        0x73, 0x2c, 0x28, 0x48, 0xcd, 0x4b, 0x51, 0x30, 0xe4, 0x32,
        0x37, 0x33, 0x31, 0xe5, 0x02, 0x00, 0xf8, 0x06, 0xaf, 0x8d,
        0x0e, 0x00, 0x00, 0x00,
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
        0x73, 0x2c, 0x28, 0x48, 0xcd, 0x4b, 0x51, 0x30, 0xe2, 0x32,
        0x33, 0x37, 0x35, 0xe1, 0x02, 0x00, 0xc6, 0x09, 0x42, 0x21,
        0x0e, 0x00, 0x00, 0x00,
    };
    CHECK(content == expected);

    // Decompress and compare
    TextFile file(filename, File::READ, File::GZIP);
    CHECK(file.readline() == "Append 1");
    CHECK(file.readline() == "7645");
    CHECK(file.readline() == "Append 2");
    CHECK(file.readline() == "6754");
    CHECK(file.readline() == "");
    CHECK(file.eof());
}

TEST_CASE("In-memory decompression") {
    auto content = std::vector<uint8_t>{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x0b, 0x49,
        0x2d, 0x2e, 0xe1, 0x32, 0x35, 0x31, 0x33, 0xe7, 0x02, 0x00, 0x8a, 0x43,
        0x5e, 0x98, 0x0a, 0x00, 0x00, 0x00,
    };

    auto decompressed = decompress_gz(reinterpret_cast<const char*>(content.data()), content.size());
    CHECK(std::string(decompressed.data(), decompressed.size()) == "Test\n5467\n");

    content[23] = 0x00;
    CHECK_THROWS_WITH(
        decompress_gz(reinterpret_cast<const char*>(content.data()), content.size()),
        "error inflating gziped memory: incorrect data check"
    );

    content[0] = 0x00;
    CHECK_THROWS_WITH(
        decompress_gz(reinterpret_cast<const char*>(content.data()), content.size()),
        "error inflating gziped memory: incorrect header check"
    );
}
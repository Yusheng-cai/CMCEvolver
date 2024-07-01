#pragma once

#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <fstream>

namespace FileSystem
{
    std::string getCurrentPath(); 

    std::string joinPath(const std::string& path1, const std::string& path2);

    std::string FileNameFromPath(const std::string& path);

    bool FileExist(const std::string& path);
};
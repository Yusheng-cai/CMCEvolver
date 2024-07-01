#pragma once 
#include "Assert.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <functional>
#include <array>
#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"

namespace StringTools
{
    using Real = CommonTypes::Real;

    // append an index to a name e.g. a.out -> a_1.out
    template <typename T>
    std::string AppendIndexToFileName(std::string str, T num, std::string delimiter=".");

    template <typename T>
    bool StringToType(std::string str, T& num);

    template <typename T>
    T StringToType(std::string str);

    template <typename T>
    std::string TypeToString(T num);

    void to_lower(std::string& str);

    // Check if a string can be converted into a numeric number (float)
    bool isNumber(std::string str);

    template<typename T>
    bool VectorStringTransform(std::vector<std::string> vecstr, std::vector<T>& output);

    void RemoveBlankInString(std::string& str);

    // Function that reads tabulated data by specified column
    template <typename T>
    void ReadTabulatedData(std::string filename, int col, std::vector<T>& data);

    // Function that reads tabulated data 
    template <typename T>
    void ReadTabulatedData(std::string filename, std::vector<std::vector<T>>& data);

    // function that reads the file extension 
    std::string ReadFileExtension(std::string filename);

    // function that reads the file name ignoring extension
    std::string ReadFileName(std::string filename, std::string delimiter=".");

    template <typename T> 
    void WriteTabulatedData(std::string filename, const std::vector<T>& data);

    template <typename T>
    void WriteTabulatedData(std::string filename, const std::vector<T>& data, const std::vector<int>& indices);
}


class ParameterPack
{
    public:
        ParameterPack():packname_("default"){};

        // Instantiate the parameter pack by name 
        ParameterPack(std::string packname):packname_(packname){};
        ~ParameterPack(){};

        enum KeyType{
            Required,
            Optional
        };

        // getters
        std::string get_packname() const {return packname_;};

        // insert into the parameterpack object
        std::string& insert(const std::string& key, const std::string& value);
        std::vector<std::string>& insert(const std::string& key, const std::vector<std::string>& value);
        ParameterPack& insert(const std::string& key, const ParameterPack& parampack);

        // Find all instances of values in Parameter Pack that matches with key
        // const function can be called by const + nonconst
        // non-const function can only be called by non-const 
        std::vector<const std::string*> findValues(const std::string& key, const KeyType) const;
        const std::string* findValue(const std::string& key, const KeyType) const;
        std::vector<const std::vector<std::string>*> findVectors(const std::string& key, const KeyType) const;
        const std::vector<std::string>* findVector(const std::string& key, const KeyType) const;
        std::vector<const ParameterPack*> findParamPacks(const std::string& key, const KeyType) const;
        const ParameterPack* findParamPack(const std::string& key, const KeyType) const;

        template <typename T>
        bool ReadNumber(const std::string& key, const KeyType, T& val) const;
        bool ReadString(const std::string& key, const KeyType, std::string& str) const;
        bool Readbool(const std::string& key, const KeyType, bool& boolean) const;

        template <typename T>
        bool ReadVectorNumber(const std::string& key, const KeyType, std::vector<T>& vecval) const;
        template <typename T>
        bool ReadVectorVectorNumber(const std::string& key, const KeyType, std::vector<std::vector<T>>& vecvecval) const;
        template <typename T, std::size_t dim>
        bool ReadVectorArrayNumber(const std::string& key, const KeyType, std::vector<std::array<T,dim>>& vecarrayval) const;
        bool ReadVectorString(const std::string& key, const KeyType, std::vector<std::string>& vecstr) const;

        template<typename T, std::size_t dim>
        bool ReadArrayNumber(const std::string& key, const KeyType, std::array<T,dim>& arrval) const;

        void print();
 
    private:
        std::multimap<std::string, std::string> value_;
        std::multimap<std::string, std::vector<std::string>> vectors_;
        std::multimap<std::string, ParameterPack> parampacks_;
        std::string packname_;
};

class TokenStream
{
    public:
        enum Status
        {
            Success, 
            Close_Brace, // }
            Open_Brace, // {
            Open_bracket, // [
            Close_bracket, // ]
            Failure,
            EndOfFile
        };

        TokenStream(std::ifstream& ifstream):ifstream_(ifstream){};

        // This is the most general way of reading tokens, we should read tokens 1 by 1 instead of reading
        // an entire line in first (former is more general) 
        Status ReadNextToken(std::string& token);
    private:   
        std::ifstream& ifstream_;
        std::stringstream line_stream_;

        std::string comment_str_="#";
};

class InputParser
{
    public:
        InputParser(){};
        ~InputParser(){};
        void ParseFile(const std::string& filename, ParameterPack& parampack);
        TokenStream::Status ParseNextToken(TokenStream& toks, ParameterPack& parampack);
        TokenStream::Status ParseParamPack(TokenStream& toks, ParameterPack& parampack);
        TokenStream::Status ParseVector(TokenStream& toks, std::vector<std::string>& vecvals);
};


template <typename T>
std::string StringTools::AppendIndexToFileName(std::string str, T num, std::string delimiter)
{
    // find out the filename 
    std::string fname = ReadFileName(str, delimiter);
    std::string extension = ReadFileExtension(str);
    std::string num_str = TypeToString(num);

    std::string ret = fname + "_" + num_str + delimiter + extension;

    return ret;
}


template <typename T>
std::string StringTools::TypeToString(T num){
    std::stringstream ss;
    ss << num;

    std::string ret;
    ss >> ret;

    return ret;
}

template <typename T>
bool StringTools::StringToType(std::string str, T& num)
{
    std::stringstream ss(str);
    ss >> num;

    return ss.fail();
}

template <typename T>
T StringTools::StringToType(std::string str)
{
    std::stringstream ss(str);

    T num;
    ss >> num;

    ASSERT((! ss.fail()), "The conversion cannot be performed with string " << str);

    return num;
}

template <typename T>
bool StringTools::VectorStringTransform(std::vector<std::string> vecstr, std::vector<T>& output)
{
    output.clear();
    for (auto str: vecstr)
    {
        T num;
        bool fail = StringTools::StringToType<T>(str, num);
        ASSERT((fail == false), "The read operation failed.");
        output.push_back(num);
    }

    return true;
}

template <typename T>
bool ParameterPack::ReadNumber(const std::string& key, const ParameterPack::KeyType keytype, T& val) const
{
    auto str = findValue(key, keytype);

    if (str != nullptr)
    {
        bool fail = StringTools::StringToType<T>(*str,val);

        ASSERT((fail == false), "The read operation failed.");

        return true;
    }

    return false;
}

template <typename T>
bool ParameterPack::ReadVectorNumber(const std::string& key, const ParameterPack::KeyType keytype, std::vector<T>& vecval) const
{
    auto vecstr = findVector(key, keytype);
    vecval.clear();

    if (vecstr != nullptr)
    {
        for (int i =0; i< vecstr->size();i++)
        {
            T val;
            bool fail = StringTools::StringToType<T>(vecstr->at(i), val);
            ASSERT((fail == false), "The read operation failed.");
            vecval.push_back(val);
        }

        return true;
    }

        
    return false;
}

template <typename T>
bool ParameterPack::ReadVectorVectorNumber(const std::string& key, const ParameterPack::KeyType keytype, std::vector<std::vector<T>>& vecvecval) const 
{
    auto vecvecstr = findVectors(key,keytype);  
    if (keytype == ParameterPack::KeyType::Required)
    {
        ASSERT((vecvecstr.size() != 0), "The required key " << key << " is not provided.");
    }

    vecvecval.clear();
    vecvecval.resize(vecvecstr.size());

    for (int i=0;i<vecvecval.size();i++)
    {
        auto& vecstr = vecvecstr[i];
        vecvecval[i].resize(vecstr->size());
        for (int j=0;j<vecstr->size();j++)
        {
            T val;
            vecvecval[i][j] = StringTools::StringToType<T>(vecstr->at(j));
        }
    }

    if (vecvecval.size() != 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <typename T, std::size_t dim>
bool ParameterPack::ReadVectorArrayNumber(const std::string& key, const ParameterPack::KeyType keytype, std::vector<std::array<T,dim>>& vecarrayval) const 
{
    std::vector<std::vector<T>> vecvecval;
    bool read = ReadVectorVectorNumber<T>(key, keytype, vecvecval);

    if ( keytype == ParameterPack::KeyType::Required)
    {
        ASSERT((read == true), "The required key " << key << " is not found.");
    }

    if (read)
    {
        vecarrayval.clear();
        vecarrayval.resize(vecvecval.size());
        for (int i=0;i<vecvecval.size();i++)
        {
            ASSERT((vecvecval[i].size() == dim), "The dimension of vector read does not match the required size of " << dim);
            for (int j=0;j<dim;j++)
            {
                vecarrayval[i][j] = vecvecval[i][j];
            }
        }

        return true;
    }
    else
    {
        return false;
    }

}

template <typename T, std::size_t dim>
bool ParameterPack::ReadArrayNumber(const std::string& key, const ParameterPack::KeyType keytype, std::array<T,dim>& arrval) const
{
    std::vector<T> vecval;
    bool vecNum = ReadVectorNumber<T>(key,keytype, vecval);


    if (vecNum == true)
    {
        ASSERT((vecval.size() == dim), "In readArrayNumber, the read vector for key= " << key << " is " << vecval.size() << " while the required size is " << dim);
        for (int i=0;i<dim;i++)
        {
            arrval[i] = vecval[i];
        }

        return true;
    }
    else
    {
        return false;
    }
}


template <typename T>
void StringTools::ReadTabulatedData(std::string filename, int col, std::vector<T>& data)
{
    std::vector<std::vector<T>> total_data;
    ReadTabulatedData(filename, total_data);

    data.clear();
    data.resize(total_data.size());
    for (int i=0;i<total_data.size();i++)
    {
        data[i] = total_data[i][col];
    }
}

template <typename T>
void StringTools::WriteTabulatedData(std::string filename, const std::vector<T>& data){
    std::ofstream ofs;
    ofs.open(filename);

    for (int i=0;i<data.size();i++){
        ofs << data[i] << "\n";
    }

    ofs.close();
}

template <typename T>
void StringTools::WriteTabulatedData(std::string filename, const std::vector<T>& data, const std::vector<int>& indices){
    ASSERT((data.size() == indices.size()), "Indices size do not match that of data.");
    std::ofstream ofs;
    ofs.open(filename);

    for (int i=0;i<data.size();i++){
        ofs << indices[i] << " " << data[i] << "\n";
    }

    ofs.close();
}


template <typename T>
void StringTools::ReadTabulatedData(std::string filename, std::vector<std::vector<T>>& data)
{
    data.clear();

    // read the filename 
    std::ifstream ifs;
    ifs.open(filename);

    ASSERT((ifs.is_open()), "The file with name " << filename << " is not opened.");

    std::string sentence;
    while(std::getline(ifs, sentence))
    {
        if (sentence.find("#") == std::string::npos)
        {
            std::stringstream ss;
            ss.str(sentence);

            std::vector<T> vecNum;
            T number;
            while (ss >> number)
            {
                vecNum.push_back(number);
            }

            data.push_back(vecNum);
        }
    }
}
#include <iostream>
#include <fstream>
#include <vector>
#include "mp_solver.h"
using namespace std;

#ifndef CSVREADER
#define CSVREADER
// This is a helper class that reads the csv for the input
// Note: while this class can handle the extra newlines in
//      end of files, it assumes that the format of the file
//      is consistent with what we have in the sample input.
//      Due to time limitation, I did implement any exception
//      handling functionality for that yet.
class csvReader
{
public:
    csvReader(string);
    ~csvReader();
    void readFile();             // Read through the target file. For each line l in the file,
                                 // store l into the vector lines.
    void readData(int);          // Read through lines and get data into expected types
    mp_config_t getData() const; // get the data in the expected type

private:
    ifstream *reader = nullptr;
    vector<string> lines;
    bool gotData;
    mp_config_t m_data;
};
#endif

csvReader::csvReader(string fileName)
    : reader{new ifstream(fileName)}, gotData{false}
{
}

csvReader::~csvReader()
{
    delete reader;
    reader = nullptr;
}

void csvReader::readFile()
{
    const int BUF_SIZE = 1024;
    char tmp[BUF_SIZE];
    while (!reader->eof())
    {
        reader->getline(tmp, BUF_SIZE);
        string line(tmp);
        if (line != "")
            lines.push_back(line);
    }
}

void csvReader::readData(int path_num = MAX_NOF_PATHS)
{
    if (lines.empty())
        this->readFile();

    m_data.nof_paths = path_num; // hardcode for now, can modify later
    // m_data.nof_pilots = MAX_NOF_PILOTS; // hardcode for now, can modify later

    int idx = 0;

    //vector<int> m, n;
    //vector<_Complex double> h;
    for (vector<string>::iterator i = lines.begin(); i < lines.end(); i++)
    {
        string tmp = *i;
        int stop_point = tmp.find(",");
        //m.push_back(stod(tmp.substr(0, stop_point)));
        m_data.m[idx] = stod(tmp.substr(0, stop_point));

        tmp = tmp.substr(stop_point + 1);
        stop_point = tmp.find(",");
        //n.push_back(stod(tmp.substr(0, stop_point)));
        m_data.n[idx] = stod(tmp.substr(0, stop_point));

        tmp = tmp.substr(stop_point + 1);
        stop_point = tmp.find("+");
        double real = stod(tmp.substr(0, stop_point));

        tmp = tmp.substr(stop_point + 1);
        stop_point = tmp.find("j");
        double imag = stod(tmp.substr(0, stop_point));
        _Complex double x = real + imag * I;
        //h.push_back(x);
        m_data.y[idx] = x;

        idx++;
    }

    m_data.nof_pilots = idx;
    gotData = true;
}

mp_config_t csvReader::getData() const
{
    if (!gotData)
    {
        cout << "didn't read data yet" << endl;
    }

    return m_data;
}
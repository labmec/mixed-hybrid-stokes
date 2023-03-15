#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

using json = nlohmann::json;

class ProblemData{
private:
    std::string fMeshName = "none";
    int fpOrder = -1;
    std::vector< std::tuple<std::string,int,double,int> > fbcvec;

public:
    ProblemData();
    void ReadJson(std::string jsonfile);
    void PrintData();
};

#endif

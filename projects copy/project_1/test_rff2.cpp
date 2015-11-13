#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;


int main() {
  const char* filename = "data_10.txt";
  ifstream inFile(filename);

  // Make sure the file stream is good
  if(!inFile) {
    cout << endl << "Failed to open file " << filename;
    return 1;
  }

  vector<double> xValues;
  vector<double> numericalValues;
  vector<double> analyticalValues;

  double x,numerical,analytical;
  while (inFile >> x >> numerical >> analytical)
  {
    xValues.push_back(x);
    numericalValues.push_back(numerical);
    analyticalValues.push_back(analytical);
  }

  return 0;
}
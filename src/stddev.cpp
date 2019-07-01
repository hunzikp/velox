#include <vector>
#include <cmath>

using namespace std;

double stddev(vector<double> scores) {
  
  double stddev;
  double sum = 0.0, mean, var = 0.0;
  
  size_t size = scores.size();
  
  
  for (int i = 0; i < size; i++){
    sum += scores[i];
  }
  
  mean = sum / size;
  
  for (int i = 0; i < size; i++){
    var += pow(scores[i] - mean, 2);
  }
  
  stddev = sqrt(var / size);
  
  return stddev;
}

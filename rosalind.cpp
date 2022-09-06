#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>

using namespace std;

void SUM3() {

  // Open file
  ifstream ifs ("3sum.txt");

  // Read number of arrays and array size
  int numArrays, arrSize;
  ifs >> numArrays >> arrSize;

  // Initialize vector of arrays
  vector<vector<int>> arrays(numArrays);
  for(int i = 0; i < numArrays; i++)
    arrays[i] = vector<int>(arrSize);

  // Read array data from file
  for(int i = 0; i < numArrays; i++)
    for(int j = 0; j < arrSize; j++)
      ifs >> arrays[i][j];

  // Close file
  ifs.close();

  for(int i = 0; i < numArrays; i++) {

    // Save indices
    map<int, vector<int> > indices;
    for(int j = 0; j < arrSize; j++)
      indices[arrays[i][j]].push_back(j);

    // Sort current array
    sort(arrays[i].begin(),arrays[i].end());

    bool brk = false;
    for(int j = 0; j < arrSize - 2; j++) {

      // Start from left and right
      int k = j + 1; int l = arrSize - 1;

      while (k != l) {

        // Sum
        int s = arrays[i][j] + arrays[i][k] + arrays[i][l];

        if (0 == s) {

          // Read and delete indices

          int i1 = indices[arrays[i][j]].back();
          indices[arrays[i][j]].pop_back();

          int i2 = indices[arrays[i][k]].back();
          indices[arrays[i][k]].pop_back();

          int i3 = indices[arrays[i][l]].back();
          indices[arrays[i][l]].pop_back();

          // Sort indices
          vector<int> ind = {i1, i2, i3};
          sort(ind.begin(),ind.end());

          // Display indices
          for (int m = 0; m < 3; m++)
            cout << ind[m] + 1 << " ";
          cout << "\n";

          brk = true; break;

        } else if (0 < s)

          // Increment right
          l--;

        else

          // Increment left
          k++;

      }

      if (brk) break;

    }

    // Indices not found for this array
    if (!brk) cout << -1 << "\n";

  }

}

int main() {
  SUM3();
}

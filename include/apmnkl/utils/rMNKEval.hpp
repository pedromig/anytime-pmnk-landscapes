/*
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; version 3
    of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Contact: http://mocobench.sourceforge.net

Authors:
    Arnaud Liefooghe <arnaud.liefooghe@lifl.fr>
    Sebastien Verel  <sebastien.verel@inria.fr>
*/

#ifndef __rMNKEval
#define __rMNKEval

/*
 *  Fitness function of the rMNK-landscapes
 *  reading a rhoMNK-landscapes instance file
 *  in c++ style
 *
 *  rhoMNK-landscapes instances can be generated with the rmnkGenerator.R
 *
 *  More information on rhoMNK-landscapes, see original paper:
 *  Verel S., Liefooghe A., Jourdan L., Dhaenens C. "Analyzing the Effect of Objective Correlation
 * on the Efficient Set of MNK-Landscapes", In Proceedings of Learning and Intelligent OptimizatioN
 * Conference (LION 5), LNCS, p. , 2011.
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace apmnkl {

namespace priv {

/// rmnk_landscapes instance evaluator
class RMNKEval {
 public:
  /*
   *  Constructor
   *
   * @param _fileName file name instance of the rho MNK-landscapes
   */
  RMNKEval(const char *_fileName) {
    load(_fileName);
  }

  /*
   *  Destructor
   */
  virtual ~RMNKEval() {
    if (links != NULL) {
      for (unsigned int n = 0; n < M; n++) {
        for (unsigned int i = 0; i < N; i++)
          delete[](links[n][i]);
        delete[] links[n];
      }

      delete[] links;
      links = NULL;
    }

    if (tables != NULL) {
      for (unsigned int n = 0; n < M; n++) {
        for (unsigned int i = 0; i < N; i++)
          delete[](tables[n][i]);
        delete[] tables[n];
      }

      delete tables;
      tables = NULL;
    }
  }

  /*
   * Compute the fitness function
   *
   * @param _solution the solution to evaluate
   * @param _objVec   the objective vector of the corresponding solution
   */
  void eval(std::vector<bool> &_solution, std::vector<double> &_objVec) {
    _objVec.resize(M);

    for (unsigned n = 0; n < M; n++)
      _objVec[n] = evalNK(n, _solution);
  }

  /*
   * to get objective space dimension
   *
   * @return dimension of the objective space
   */
  unsigned getM() {
    return M;
  }

  /*
   * to get bitstring size
   *
   * @return dimension of the bitstring
   */
  unsigned getN() {
    return N;
  }

  /*
   * to get epistasis degree (K)
   *
   * @return epistasis degree K
   */
  unsigned getK() {
    return K;
  }

  /*
   * to get the correlation between each tuple of contributions
   *
   * @return parameter rho
   */
  double getRho() {
    return rho;
  }

 protected:
  // correlation between contributions
  double rho;

  // number of objective functions
  unsigned M;

  // size of the bit string
  unsigned N;

  // number of interactions between variables (epistasis)
  unsigned K;

  // the M tables of contributions
  double ***tables;

  // the M links description
  unsigned ***links;

  /***********************************************
   *
   * Load the file of a rMNK-landscapes instance
   *
   * @param fileName file name instance of the rMNK-landscapes
   *
   ***********************************************/
  virtual void load(const char *_fileName) {
    std::fstream file;
    file.open(_fileName, std::ios::in);

    if (file.is_open()) {
      std::string s;

      // read the commentaries
      std::string line;

      file >> s;
      while (s[0] == 'c') {
        getline(file, line, '\n');
        file >> s;
      }

      // read the parameters
      if (s.compare("p") != 0)
        std::cerr << "Error RMNKEval.load: expected line beging by \"p\" at " + s + " in "
                  << _fileName;

      file >> s;

      if (s.compare("rMNK") != 0)
        std::cerr << "Error RMNKEval.load: type rMNK expected at " + s + " in " << _fileName;

      // effective read of the parameters
      file >> rho >> M >> N >> K;

      init();

      file >> s;
      // read the links
      if (s.compare("p") != 0)
        std::cerr << "Error RMNKEval.load: expected line beging by \"p\" at " + s + " in "
                  << _fileName;

      file >> s;
      if (s.compare("links") == 0)
        loadLinks(file);
      else
        std::cerr << "Error RMNKEval.load: line with \"links\" expected at " + s + " in "
                  << _fileName;

      // read the tables of contributions
      file >> s;
      if (s.compare("p") != 0)
        std::cerr << "Error RMNKEval.load: expected line beging by \"p\" at " + s + " in "
                  << _fileName;

      file >> s;

      if (s.compare("tables") == 0)
        loadTables(file);
      else
        std::cerr << "Error RMNKEval.load: line with \"tables\" expected at " + s + " in "
                  << _fileName;

      file.close();
    } else
      std::cerr << "Error RMNKEval.load: impossible to open file " << _fileName;
  };

  /***********************************************
   *
   * Initialization of the different tables and epistasis links
   *
   ***********************************************/
  void init() {
    links = new unsigned **[M];
    tables = new double **[M];

    for (unsigned n = 0; n < M; n++) {
      links[n] = new unsigned *[N];
      tables[n] = new double *[N];

      for (unsigned i = 0; i < N; i++) {
        tables[n][i] = new double[1 << (K + 1)];
        links[n][i] = new unsigned[K + 1];
      }
    }
  }

  /***********************************************
   *
   * Load the epistasis links from file
   *
   * @param _file open file of the instance
   *
   ***********************************************/
  void loadLinks(std::fstream &_file) {
    unsigned n, i, j;

    for (i = 0; i < N; i++)
      for (j = 0; j < K + 1; j++)
        for (n = 0; n < M; n++)
          _file >> links[n][i][j];
  }

  /***********************************************
   *
   * Load the tables of contribution
   *
   * @param _file open file of the instance
   *
   ***********************************************/
  void loadTables(std::fstream &_file) {
    unsigned n, i;
    int j;

    for (i = 0; i < N; i++)
      for (j = 0; j < (1 << (K + 1)); j++)
        for (n = 0; n < M; n++)
          _file >> tables[n][i][j];
  }

  /***********************************************
   *
   * Fitness function of a single-objective NK-landscapes
   *
   * @param _numObj the objective fuction to consider
   * @param _sol the solution to evaluate
   *
   ***********************************************/
  double evalNK(unsigned _numObj, std::vector<bool> &_sol) {
    double accu = 0.0;

    for (unsigned int i = 0; i < N; i++)
      accu += tables[_numObj][i][sigma(_numObj, _sol, static_cast<int>(i))];

    return accu / static_cast<double>(N);
  }

  /***********************************************
   *
   * Extract epistatic links of the fitness contribution i
   *
   * @param _numObj the objective function to consider
   * @param _sol the solution to evaluate
   * @param _i bit of the contribution
   *
   * **********************************************/
  unsigned int sigma(unsigned _numObj, std::vector<bool> &_sol, int _i) {
    unsigned int n = 1;
    unsigned int accu = 0;

    for (unsigned int j = 0; j < K + 1; j++) {
      if (_sol[links[_numObj][_i][j]] == 1)
        accu = accu | n;

      n = n << 1;
    }

    return accu;
  }
};

}  // namespace priv
}  // namespace pmnk
#endif  // __rMNKEval

#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace Rcpp;



// check the last ending characters in a word
bool hasEnding (std::string const &fullString, std::string const &ending) {
  if (fullString.length() >= ending.length()) {
    return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
  } else {
    return false;
  }
}


//' read ms-formatted MaCS output into a matrix of haplotypes.
//'
//' This is implemented in Rcpp for speed.  It is designed to read a single replicate coalescent
//' simulation in macs' msformatter format and convert the 0/1's to 1/2's and return in a big integer
//' matrix of haplotypes that can then be put together as desired later. (i.e. different pops and into
//' different diploid individuals)
//' @param Input the path to the mac | msformatter output file to read in. Note that
//' tilde expansion will not occur on this path
//' @export
//' @keywords internal
// [[Rcpp::export]]
List read_macs_output(CharacterVector Input) {
  int i,j;
  int N, S;
  std::string tempstr;
  std::string line;
  std::string word;
  std::string fname = Rcpp::as<std::string>(Input);



  // open the file to read
  std::ifstream infile (fname);

  // first read the first line and parse it to get the second token which gives us the number of haplotypes
  std::getline(infile, line);
  std::stringstream ss(line);
  ss >> word;
  if(!hasEnding(word, "macs")) {
    stop("macs output file does not seem to start with the right command line on first line");
  }
  ss >> word;
  N = std::stoi(word);  // number of haplotypes/sequences that we will have here

  // now read lines down to the segsites
  for(j=0;j<4;j++) std::getline(infile, line);
  std::stringstream ss2(line);
  ss2 >> word;
  Rcpp::Rcout << word << std::endl;
  if(word != "segsites:") {
    stop("Sorry, I expected to find the segsites line here....");
  }
  ss2 >> word;
  S = std::stoi(word);

  // now read the line that has the positions:
  std::getline(infile, line);
  std::stringstream ss3(line);
  ss3 >> word;
  if(word != "positions:") {
    stop("Expected to find the positions line here...");
  }
  NumericVector Pos(S);  // declare space to return this
  for(j=0;j<S;j++) {
    ss3 >> word;
    Pos[j] = std::stod(word);
  }

  // and finally we read in the haplotypes and put them in a matrix with seqs in rows and segregating sites in columns
  IntegerMatrix Mat(N, S);
  for(j=0;j<N;j++) {
    std::getline(infile, line);
    for(i=0;i<S;i++) {
      Mat(j,i) = line.at(i) - '0' + 1;
    }
  }

  Rcpp::Rcout << "Num Segsites: " << S << "  Num Seqs: " << N << std::endl;
  List ret;
  ret["pos"] = Pos;
  ret["mat"] = Mat;
  return(ret);
}

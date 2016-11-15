#include "CRprintInfo.h"

void CRprintInfo(TString val, TString text) {
  TString line;
  for(unsigned int i=0; i<88; i++) line += "-";
  std::cout << boost::format("%-15s # %-70s") % val % text << std::endl;
  std::cout << line << std::endl;
}

void CRprintInfo(bool val, TString text) {
  TString val_str;
  if(val) {
    val_str = "yes|0==   |no";
  }else {
    val_str = "yes|   ==0|no";
  }
  TString line;
  for(unsigned int i=0; i<90; i++) line += "-";
  std::cout << boost::format("%-15s --------  %-70s") % val_str % text << std::endl;
  //std::cout << line << std::endl;
}

void CRprintInfo(float val, TString text) {
  TString line;
  for(unsigned int i=0; i<90; i++) line += "-";
  std::cout << boost::format("%20.4f ---  %-70s") % val % text << std::endl;
  //std::cout << line << std::endl;
}

void CRprintInfo(int val, TString text) {
  TString line;
  for(unsigned int i=0; i<90; i++) line += "-";
  std::cout << boost::format("%20d ---  %-70s") % val % text << std::endl;
  //std::cout << line << std::endl;
}

void CRprintInfo(unsigned int val, TString text) {
  TString line;
  for(unsigned int i=0; i<90; i++) line += "-";
  std::cout << boost::format("%20d ---  %-70s") % val % text << std::endl;
  //std::cout << line << std::endl;
}

#pragma once

#include<cstdlib>
#include<cstdio>
#include<string>

#define RED_STRING(string) "\x1b[31m" string "\x1b[0m"
#define GREEN_STRING(string) "\x1b[32m" string "\x1b[0m"

void pass(const char *filename, const char *funcname, const std::string &comment) {
  std::string str_file(filename), str_func(funcname);
  std::string out = str_file + "\t" + str_func;
  if(comment.size() > 0){
    out + " : ";
    out += comment;
  }
  out += " \t PASS\n";
  fprintf(stderr, GREEN_STRING("%s"), out.c_str());
  return;
}

void pass(const char *filename, const char *funcname, const char *comment){
  std::string comment_str(comment);
  pass(filename, funcname, comment_str);
  return;
}

void pass(const char *filename, const char *funcname){
  std::string tmp;
  pass(filename, funcname, tmp);
  return;
}

void failed(const char *filename, const char *funcname, const std::string &comment) {
  std::string str_file(filename), str_func(funcname);
  std::string out = str_file + "\t" + str_func;
  if(comment.size() > 0){
    out + " : ";
    out += comment;
  }
  out += " \t FAILED\n";
  fprintf(stderr, RED_STRING("%s"), out.c_str());
  return;
}

void failed(const char *filename, const char *funcname, const char *comment){
  std::string comment_str(comment);
  failed(filename, funcname, comment_str);
  return;
}

void failed(const char *filename, const char *funcname){
  std::string tmp;
  failed(filename, funcname, tmp);
  return;
}

template<typename TYPE>
std::string get_type(void) { return typeid(TYPE).name(); }
template <> std::string get_type<double>(void) { return "double"; }
template <> std::string get_type<float>(void) { return "float"; }

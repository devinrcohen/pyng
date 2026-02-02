//
// Created by devinrcohen on 1/29/26.
//

#include "montecarlo.hpp"
#include <unordered_map>
#include <random>

using namespace std;

RandomComponent::RandomComponent(std::string compname, float nom, float tol, DistributionType type) :
nom(nom), tol(tol), dist(type), name(compname){}

void test() {
    unordered_map<int, int> m;
}
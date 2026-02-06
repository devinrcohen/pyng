//
// Created by devinrcohen on 1/29/26.
//

#include "montecarlo.hpp"
#include <unordered_map>
#include <random>

using namespace std;

RandomComponent::RandomComponent(std::string compname, float nom, float tol, DistributionType type) :
name(compname), nom(nom), tol(tol), dist(type) {}

void test() {
    unordered_map<int, int> m;
}
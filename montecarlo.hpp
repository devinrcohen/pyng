//
// Created by devinrcohen on 1/29/26.
//

#pragma once
#include <string>
#include <unordered_map>

enum DistributionType {
    Uniform,
    Gaussian,
    Normal
};

class RandomComponent {
public:
    RandomComponent(std::string, float, float, DistributionType);
    float Nom() const noexcept { return nom; }
    float Tol() const noexcept { return tol; }
    DistributionType Dist() const noexcept { return dist; }
    std::string Name() const noexcept { return name; }
private:
    std::string name;
    DistributionType dist;
    float nom; // std. dev for normal/Gaussian
    float tol;
};
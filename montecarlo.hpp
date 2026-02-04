//
// Created by devinrcohen on 1/29/26.
//

#pragma once
#include <string>
#include <unordered_map>
#include <random>
#include <tuple>

enum DistributionType {
    Uniform,
    Gaussian,
    Normal
};

class RandomComponent {
public:
    RandomComponent(std::string, float, float, DistributionType);
    [[nodiscard]] float Nom() const noexcept { return nom; }
    [[nodiscard]] float Tol() const noexcept { return tol; }
    [[nodiscard]] DistributionType Dist() const noexcept { return dist; }
    [[nodiscard]] std::string Name() const noexcept { return name; }
    static void component(std::tuple<std::string, std::string, std::string, std::string>);
private:
    std::string name;
    DistributionType dist;
    float nom; // std. dev for normal/Gaussian
    float tol;
};
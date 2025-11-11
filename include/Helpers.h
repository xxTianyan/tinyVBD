//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef HELPER_H
#define HELPER_H

#include <iostream>

template<typename T>
void print_vector(const std::vector<T>& vec) {
    for (const auto& elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}


#endif //HELPER_H

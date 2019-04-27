#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdio>
#include <set>
#include <cassert>

const long long MATRIX_SIZE_MAX = 100000000;

class BackpackSolution {
public:
    BackpackSolution() {
        Weight = 0;
        Value = 0;
    }
    
    BackpackSolution(const BackpackSolution&& rvalue) {
        this->Value = rvalue.Value;
        this->Weight = rvalue.Weight;
        this->TakenItems = std::move(rvalue.TakenItems);
    }

    void Merge(const BackpackSolution& other) {
        this->Weight += other.Weight;
        this->Value += other.Value;
        for (size_t i = 0; i < other.TakenItems.size(); ++i) {
            this->TakenItems.push_back(other.TakenItems[i]);
        }
    }

    int Value;
    int Weight;
    std::vector<int> TakenItems;

private:
    // копировать запрещаес
    BackpackSolution(const BackpackSolution&) = delete;
};

struct Item {
    int Weight = 0;
    int Value = 0;
    float Usability = 0.;
    int Index = 0;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
// Вспомогательные функции
/////////////////////////////////////////////////////////////////////////////////////////////////

// Главная функция для поиска решения по матрице "полезности"
static void find_solution_recursive(int i, int w,
        const std::vector<std::vector<int>>& solutionMatrix,
        const std::vector<Item>& items,
        BackpackSolution& solution) {
    if (solutionMatrix[i][w] == 0) {
        return;
    }
    if (solutionMatrix[i - 1][w] == solutionMatrix[i][w]) {
        find_solution_recursive(i - 1, w, solutionMatrix, items, solution);
    } else {
        find_solution_recursive(i - 1, w - items[i - 1].Weight, solutionMatrix, items, solution);
        solution.TakenItems.push_back(items[i - 1].Index);
        solution.Weight += items[i - 1].Weight;
    }
    
}

// вызов главной функции без кучи ненужных аргументов
static BackpackSolution find_solution(const std::vector<std::vector<int>>& solutionMatrix,
        const std::vector<Item>& items) {
    BackpackSolution solution;
    find_solution_recursive(solutionMatrix.size() - 1, solutionMatrix[0].size() - 1,
            solutionMatrix, items, solution);
    return solution;
}

// функция сравнения двух Item для жадного алгоритма
static bool itemUsabilityCmp(const Item& first, const Item& second) {
    return first.Usability > second.Usability ||
        (first.Usability == second.Usability && first.Value > second.Value);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

BackpackSolution solve_backpack_greedy(std::vector<Item>& items, int maxWeight) {
    BackpackSolution solution;
    std::sort(items.begin(), items.end(), itemUsabilityCmp);
    for (size_t i = 0; i < items.size(); ++i) {
        Item& item = items[i];
        if (solution.Weight + item.Weight <= maxWeight) {
            solution.Value += item.Value;
            solution.Weight += item.Weight;
            solution.TakenItems.push_back(item.Index);
        }
    }
    return solution;
}

BackpackSolution solve_backpack_dymamic(const std::vector<Item>& items, int maxWeight) {
    std::vector<std::vector<int>> solutionMatrix(items.size() + 1,
            std::vector<int>(maxWeight + 1, 0));

    for (size_t i = 1; i < items.size() + 1; ++i) {
        for (int w = 1; w < maxWeight + 1; ++w) {
            if (w >= items[i - 1].Weight) {
                solutionMatrix[i][w] = std::max(solutionMatrix[i - 1][w], 
                        solutionMatrix[i - 1][w - items[i - 1].Weight] + items[i - 1].Value);
            } else {
                solutionMatrix[i][w] = solutionMatrix[i - 1][w];
            }
        }
    }
    BackpackSolution solution = find_solution(solutionMatrix, items);
    solution.Value = solutionMatrix.back().back();
    return solution;
}

int main(int argc, char* argv[]) {
    // Эта ересь, чтобы можно было в аргументах ком. строки передавать файл на вход и выход
    std::ifstream input;
    std::ofstream output;
    if (argc > 1) {
        input = std::ifstream(argv[1]);
        std::cin.rdbuf(input.rdbuf());
    }
    if (argc > 2) {
        output = std::ofstream(argv[2]);
        std::cout.rdbuf(output.rdbuf());
    }
    
    int itemsCount, maxWeight;
    std::cin >> itemsCount >> maxWeight;
    std::vector<Item> items(itemsCount);
    
    for (size_t i = 0; i < items.size(); ++i) {
        std::cin >> items[i].Value >> items[i].Weight;
        items[i].Usability = static_cast<float>(items[i].Value) / items[i].Weight;
        items[i].Index = i + 1;
    }

    int tail = MATRIX_SIZE_MAX / itemsCount;
    BackpackSolution backpack = solve_backpack_greedy(items, maxWeight - tail);
    
    // смотрим, сколько не хватает до полной вместимости
    int tailMaxWeight = maxWeight - backpack.Weight;
    assert(tailMaxWeight >= 0);
    std::vector<Item> tailItems;
    std::set<int> takenInds(backpack.TakenItems.begin(), backpack.TakenItems.end());
    // наберем только те, что вмещаются
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].Weight <= tailMaxWeight && takenInds.count(items[i].Index) == 0) {
            tailItems.push_back(items[i]);
        }
    }
    // соединяем 2 решения
    backpack.Merge(solve_backpack_dymamic(tailItems, tailMaxWeight));

    std::cout << backpack.Value << std::endl;
    for (size_t i = 0; i < backpack.TakenItems.size(); ++i) {
        std::cout << backpack.TakenItems[i];
        if (i != backpack.TakenItems.size() - 1) {
            std::cout << " ";
        } else {
            std::cout << std::endl;
        }
    }
    return 0;
}


#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <string>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class SingleQubit {
private:
    std::complex<double> alpha;
    std::complex<double> beta;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dis;

public:
    SingleQubit() : alpha(1.0, 0.0), beta(0.0, 0.0), gen(rd()), dis(0.0, 1.0) {
        normalizeState();
    }

    void normalizeState() {
        double norm = std::sqrt(std::norm(alpha) + std::norm(beta));
        if (norm > 1e-10) { 
            alpha /= norm;
            beta /= norm;
        }
        else {
            alpha = std::complex<double>(1.0, 0.0);
            beta = std::complex<double>(0.0, 0.0);
        }
    }
    void hadamard() {
        const double sqrt2 = std::sqrt(2.0);
        std::complex<double> new_alpha = (alpha + beta) / sqrt2;
        std::complex<double> new_beta = (alpha - beta) / sqrt2;
        alpha = new_alpha;
        beta = new_beta;
        normalizeState();
    }
    void pauliX() {
        std::swap(alpha, beta);
    }

    void pauliZ() {
        beta = -beta;
    }
    void pauliY() {
        std::complex<double> i(0.0, 1.0);
        std::complex<double> new_alpha = -i * beta;
        std::complex<double> new_beta = i * alpha;
        alpha = new_alpha;
        beta = new_beta;
        normalizeState(); 
    }

    void rotateX(double theta) {
        double cos_half = std::cos(theta / 2);
        double sin_half = std::sin(theta / 2);
        std::complex<double> new_alpha = cos_half * alpha - std::complex<double>(0, 1) * sin_half * beta;
        std::complex<double> new_beta = -std::complex<double>(0, 1) * sin_half * alpha + cos_half * beta;
        alpha = new_alpha;
        beta = new_beta;
        normalizeState();
    }

    void rotateY(double theta) {
        double cos_half = std::cos(theta / 2);
        double sin_half = std::sin(theta / 2);
        std::complex<double> new_alpha = cos_half * alpha - sin_half * beta;
        std::complex<double> new_beta = sin_half * alpha + cos_half * beta;
        alpha = new_alpha;
        beta = new_beta;
        normalizeState();
    }

    void rotateZ(double theta) {
        double cos_half = std::cos(theta / 2);
        double sin_half = std::sin(theta / 2);
        std::complex<double> phase_0(cos_half, -sin_half);
        std::complex<double> phase_1(cos_half, sin_half);

        alpha *= phase_0;
        beta *= phase_1;
        normalizeState();
    }

    void rotate(double theta, double nx, double ny, double nz) {
       
        double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (norm > 1e-10) {
            nx /= norm;
            ny /= norm;
            nz /= norm;
        }
        else {
            nx = 0; ny = 0; nz = 1;
        }
        double cos_half_theta = std::cos(theta / 2);
        double sin_half_theta = std::sin(theta / 2);
        std::complex<double> a(cos_half_theta, -sin_half_theta * nz);
        std::complex<double> b(-sin_half_theta * ny, -sin_half_theta * nx);
        std::complex<double> c(sin_half_theta * ny, -sin_half_theta * nx);
        std::complex<double> d(cos_half_theta, sin_half_theta * nz);
        std::complex<double> new_alpha = a * alpha + b * beta;
        std::complex<double> new_beta = c * alpha + d * beta;

        alpha = new_alpha;
        beta = new_beta;
        normalizeState();
    }

    int measure() {
        double prob_zero = std::norm(alpha);
        double random_val = dis(gen);

        if (random_val < prob_zero) {
            alpha = std::complex<double>(1.0, 0.0);
            beta = std::complex<double>(0.0, 0.0);
            return 0;
        }
        else {
            alpha = std::complex<double>(0.0, 0.0);
            beta = std::complex<double>(1.0, 0.0);
            return 1;
        }
    }

    std::pair<double, double> getProbabilities() const {
        return { std::norm(alpha), std::norm(beta) };
    }

    std::pair<std::complex<double>, std::complex<double>> getState() const {
        return { alpha, beta };
    }

    void setState(std::complex<double> new_alpha, std::complex<double> new_beta) {
        alpha = new_alpha;
        beta = new_beta;
        normalizeState();
    }

    void printState() const {
        std::pair<double, double> probs = getProbabilities();
        double prob0 = probs.first;
        double prob1 = probs.second;

        std::cout << "Состояние кубита: " << std::fixed << std::setprecision(4);
        if (std::abs(alpha.imag()) < 1e-10) {
            std::cout << "(" << alpha.real() << ")|0⟩ + ";
        }
        else if (alpha.imag() < 0) {
            std::cout << "(" << alpha.real() << " - " << std::abs(alpha.imag()) << "i)|0⟩ + ";
        }
        else {
            std::cout << "(" << alpha.real() << " + " << alpha.imag() << "i)|0⟩ + ";
        }

        if (std::abs(beta.imag()) < 1e-10) {
            std::cout << "(" << beta.real() << ")|1⟩" << std::endl;
        }
        else if (beta.imag() < 0) {
            std::cout << "(" << beta.real() << " - " << std::abs(beta.imag()) << "i)|1⟩" << std::endl;
        }
        else {
            std::cout << "(" << beta.real() << " + " << beta.imag() << "i)|1⟩" << std::endl;
        }

        std::cout << "Вероятности: |0⟩: " << prob0 * 100 << "%, |1⟩: " << prob1 * 100 << "%" << std::endl;
    }
};

int main() {
    setlocale(LC_ALL, "");

    std::cout << "Квантовый калькулятор на одном кубите" << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    SingleQubit qubit;
    std::string input;

    std::cout << "Начальное состояние кубита (|0⟩):" << std::endl;
    qubit.printState();

    while (true) {
        std::cout << "\nДоступные команды:" << std::endl;
        std::cout << "h     - применить гейт Адамара (H)" << std::endl;
        std::cout << "x     - применить гейт Паули-X" << std::endl;
        std::cout << "y     - применить гейт Паули-Y" << std::endl;
        std::cout << "z     - применить гейт Паули-Z" << std::endl;
        std::cout << "rx N  - повернуть на N градусов вокруг оси X" << std::endl;
        std::cout << "ry N  - повернуть на N градусов вокруг оси Y" << std::endl;
        std::cout << "rz N  - повернуть на N градусов вокруг оси Z" << std::endl;
        std::cout << "r N X Y Z - повернуть на N градусов вокруг произвольной оси [X,Y,Z]" << std::endl;
        std::cout << "m     - измерить кубит" << std::endl;
        std::cout << "state - показать текущее состояние кубита" << std::endl;
        std::cout << "set a1 a2 b1 b2 - установить состояние (a1+a2i)|0⟩ + (b1+b2i)|1⟩" << std::endl;
        std::cout << "q     - выход" << std::endl;
        std::cout << "> ";

        std::getline(std::cin >> std::ws, input);
        std::istringstream iss(input);
        std::string command;
        iss >> command;

        try {
            if (command == "q") {
                break;
            }
            else if (command == "h") {
                qubit.hadamard();
                std::cout << "Применен гейт Адамара (H)" << std::endl;
            }
            else if (command == "x") {
                qubit.pauliX();
                std::cout << "Применен гейт Паули-X" << std::endl;
            }
            else if (command == "y") {
                qubit.pauliY();
                std::cout << "Применен гейт Паули-Y" << std::endl;
            }
            else if (command == "z") {
                qubit.pauliZ();
                std::cout << "Применен гейт Паули-Z" << std::endl;
            }
            else if (command == "rx") {
                double angle;
                if (!(iss >> angle)) {
                    throw std::runtime_error("Требуется числовое значение угла");
                }
                qubit.rotateX(angle * M_PI / 180.0); 
                std::cout << "Применен поворот на " << angle << " градусов вокруг оси X" << std::endl;
            }
            else if (command == "ry") {
                double angle;
                if (!(iss >> angle)) {
                    throw std::runtime_error("Требуется числовое значение угла");
                }
                qubit.rotateY(angle * M_PI / 180.0);
                std::cout << "Применен поворот на " << angle << " градусов вокруг оси Y" << std::endl;
            }
            else if (command == "rz") {
                double angle;
                if (!(iss >> angle)) {
                    throw std::runtime_error("Требуется числовое значение угла");
                }
                qubit.rotateZ(angle * M_PI / 180.0);
                std::cout << "Применен поворот на " << angle << " градусов вокруг оси Z" << std::endl;
            }
            else if (command == "r") {
                double angle, nx, ny, nz;
                if (!(iss >> angle >> nx >> ny >> nz)) {
                    throw std::runtime_error("Требуются числовые значения: угол и координаты оси [X Y Z]");
                }
                qubit.rotate(angle * M_PI / 180.0, nx, ny, nz);
                std::cout << "Применен поворот на " << angle << " градусов вокруг оси ["
                    << nx << ", " << ny << ", " << nz << "]" << std::endl;
            }
            else if (command == "m") {
                int result = qubit.measure();
                std::cout << "Результат измерения: |" << result << "⟩" << std::endl;
            }
            else if (command == "state") {
            }
            else if (command == "set") {
                double a1, a2, b1, b2;
                if (!(iss >> a1 >> a2 >> b1 >> b2)) {
                    throw std::runtime_error("Требуются 4 числовых значения для установки состояния");
                }
                std::complex<double> new_alpha(a1, a2);
                std::complex<double> new_beta(b1, b2);
                qubit.setState(new_alpha, new_beta);
                std::cout << "Установлено новое состояние кубита" << std::endl;
            }
            else {
                std::cout << "Неизвестная команда: " << command << std::endl;
            }
        }
        catch (const std::exception& e) {
            std::cout << "Ошибка: " << e.what() << std::endl;
        }

        qubit.printState();
    }

    std::cout << "Выход из квантового калькулятора" << std::endl;
    return 0;
}

#pragma once
#include <chrono>
#include <iostream>

class TimeLog {
public:
    TimeLog(std::string const &title): _title(title) {
        restart();
    }

    void restart() {
        begin = std::chrono::steady_clock::now();
        lastSubStep = begin;
    }

    void logSubStep(std::string subTitle) {
        auto now = std::chrono::steady_clock::now();
        std::cout << "    [Time log] " << _title << " > " << subTitle << " = " << std::chrono::duration_cast<std::chrono::milliseconds>(now - lastSubStep).count() / 1000. <<"s."<< std::endl;
        lastSubStep = now;
    }

    void logTotalTime() {
        auto now = std::chrono::steady_clock::now();
        std::cout << "[Time log] " << _title << " > " << "TOTAL TIME = " << std::chrono::duration_cast<std::chrono::milliseconds>(now - begin).count() / 1000. << "s."<< std::endl;
    }
private:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point lastSubStep;

    std::string _title;
};


    #include <iostream>
    #include <vector>
    #include <algorithm>
    #include <tr1/random>

    #include <omp.h>  // For omp_get_wtime() only.

    const int ITERATIONS = 10;
    const int N = 100*1000*1000;
    const int SEED = 42;

    struct Item1
    {
        int i1, i2, i3, i4, i5;

        Item1() : i1(0), i2(0), i3(0), i4(0), i5(0)
        {}

        Item1(int i) : i1(i1), i2(i), i3(i), i4(i), i5(i)
        {}

        bool operator<(Item1 const & other) const
        {
            return i1 < other.i1;
        }
    };

    struct Item2
    {
        int i1, i2, i3, i4, i5;

        Item2() : i1(0), i2(0), i3(0), i4(0), i5(0)
        {}

        Item2(int i) : i1(i1), i2(i), i3(i), i4(i), i5(i)
        {}

        virtual
        bool operator<(Item2 const & other) const
        {
            return i1 < other.i1;
        }
    };

    int main()
    {
        double start = 0;
        double timeAvg = 0;
        double timeDev = 0;
        std::vector<double> times;

        std::cout << "Parameters\n";
        std::cout << "    # iterations: " << ITERATIONS << "\n";
        std::cout << "    # items     : " << N << "\n";
        std::cout << "    seed        : " << SEED << "\n\n";

        // Generate random input.
        std::cout << "Generating input.\n";
        start = omp_get_wtime();
        std::tr1::mt19937 eng(SEED);
        std::tr1::uniform_int<int> unif;
        std::vector<int> nums;
        nums.reserve(N);
        for (int i = 0; i < N; ++i)
            nums.push_back(unif(eng));
        std::cout << "    time " << omp_get_wtime() - start << " s\n\n";

        // Sort with inlining.
        times.clear();
        timeAvg = 0;
        timeDev = 0;
        std::cout << "std::sort with inlining\n";
        for (int i = 0; i < ITERATIONS + 1; ++i)
        {
            std::vector<Item1> items(nums.begin(), nums.end());
            start = omp_get_wtime();
            std::sort(items.begin(), items.end());
            if (i > 0)
              times.push_back(omp_get_wtime() - start);
        }
        for (unsigned i = 0; i < times.size(); ++i)
            timeAvg += times[i];
        timeAvg /= times.size();
        for (unsigned i = 0; i < times.size(); ++i)
            timeDev += (timeAvg - times[i]) * (timeAvg - times[i]);
        timeDev /= times.size();
        timeDev = sqrt(timeDev);
        std::cout << "    time avg " << timeAvg << " s dev " << timeDev << "\n\n";

        // Sorting with virtual operator<().
        times.clear();
        timeAvg = 0;
        timeDev = 0;
        std::cout << "std::sort without inlining\n";
        for (int i = 0; i < ITERATIONS + 1; ++i)
        {
            std::vector<Item2> items(nums.begin(), nums.end());
            start = omp_get_wtime();
            std::sort(items.begin(), items.end());
            if (i > 0)
              times.push_back(omp_get_wtime() - start);
        }
        for (unsigned i = 0; i < times.size(); ++i)
            timeAvg += times[i];
        timeAvg /= times.size();
        for (unsigned i = 0; i < times.size(); ++i)
            timeDev += (timeAvg - times[i]) * (timeAvg - times[i]);
        timeDev /= times.size();
        timeDev = sqrt(timeDev);
        std::cout << "    time avg " << timeAvg << " s dev " << timeDev << "\n";

        return 0;
    }

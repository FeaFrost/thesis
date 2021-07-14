
#include <iostream> // Работа с консолью
#include <fstream>  // Работа с файлами
#include <cstdlib>  // Для exit()
#include <iomanip>  // Для ровного заполнения и длины чисел

#include <random> //
#include <ctime>  // Для вихря
#include <time.h> // Для случайного распределения Пуасона
#include <string> // Строки
#include <vector> // Для дин массивов
#include "random.h"

//~----------------------------------------FUNCTIONS-----------------------------------------
long long mTwist_int(long long from, long long till)
{ // Вихрь Мерсенна для целых чисел в диапазоне from-till
    static std::mt19937 gen(time(0));
    std::uniform_int_distribution<> uid(from, till);
    return uid(gen);
}
double mTwist_real(double from, double till)
{ // Вихрь Мерсенна для действительных чисел в диапазоне from-till
    static std::mt19937 gen(time(0));
    std::uniform_real_distribution<> uid(from, till);
    return uid(gen);
}
void randomPuasonTime()
{
    time_t *td = new time_t;
    time(td);
    unsigned long init[4] = {0x123 + unsigned(*td), 0x234, 0x345, 0x456}, length = 4;
    init_by_array(init, length);
}
void stopClosing()
{
    std::cin.clear();
    std::cin.ignore(32767, '\n');
    std::cin.get();
}
long long llProbTy(std::string prStr)
{
    int leng = prStr.length() - 2;
    long long llPr = (long long)round(((std::stod(prStr)) * pow(10, leng)));
    return (llPr);
}

int** matrxDyn(int n, int m)
{
    int **A;// Создаем матрицу введенной размерности
    A = new int *[n]; // через массив указателей
    for (int i = 0; i < n; i++) {
       A[i] = new int [m];
    }
    return A;
}

template <typename Type>
std::string to_str(const Type &t)
{
    std::ostringstream os;
    os << t;
    return os.str();
}

using namespace std;
//~------------------------------------------------------------------------------------------




//~-------------------------------------------MAIN-------------------------------------------
int main() {
    int conStepsZ{1}, conSteps{100};

    string filepath2 = "results/All_steps.txt"; // путь к файлу с изменяемым именем
    ofstream outS(filepath2, std::ios::out);
    if (!outS)
    {
        std::cerr << endl
                  << "All_steps.txt couldn't be opened for writing!" << endl;
    }
    outS.close(); //Закрываем файл

    randomPuasonTime(); // Для случайного распределения Пуасона




//* ANCHOR 1 -- Читаем настройки из файла
    std::cout << "1. Importing Data from options.txt.." << endl;
    int dnaSize{1000}, conSteps2{100}, conInt{95}, nichNumb;
    string probTyS{"1"};
    double mutProbTy, jumpProbTy;

    ifstream inf("options.txt");
    if (!inf){
        std::cerr << "ERR. options.txt NOT FOUND or couldn't be opened for reading!" << endl;
        stopClosing();
        exit(1);
    }
    else{
        string blank{""};
        inf >> blank >> dnaSize >> blank >> mutProbTy >> blank >> conInt >> blank >> conSteps2 >> blank >>
            nichNumb >> blank >> jumpProbTy;

        probTyS = to_str(mutProbTy);

        inf.close(); //Закрываем файл
        std::cout << " SUCCESS." << endl;
    }

    //^ Вывод считанных опций
    std::cout << endl
              << "             Your options" << endl
              << "   ---------------------------------" << endl
              << "   DNA size                  : " << dnaSize << endl
              << "   Mutation PROBABILITY      : " << mutProbTy << endl
              << "   Confidence interval       : " << conInt << endl
              << "   Steps for c.interval      : " << conSteps2 << endl
              << "   Number of niches          : " << nichNumb << endl
              << "   Jump to new niche PROB-TY : " << jumpProbTy << endl
              << "   ---------------------------------" << endl;
//* 1 --------------------------------------------




//* ANCHOR 2 -- Создаём массив из случайных нуклеотидов
    std::cout << endl
              << "2. Initialization..." << endl;

    int d{};
    int digiotide[dnaSize]{};   
    int digiotideFour[dnaSize]{};
    
    for (int nuclExact = 0; nuclExact < dnaSize; nuclExact++)
    {
        d = mTwist_int(0, 3);
        digiotide[nuclExact] = d;
        digiotideFour[nuclExact] = 4;
    }

    std::cout << "3. DNA synthesis was successful." << endl;
//* 2 --------------------------------------------




//* ANCHOR 3 -- Начинаем распределение по нишам
    std::cout << "4. Starting JUMPING.." << endl;




//* ANCHOR 3.1 -- Вводим переменные для подсчёта занятых ниш
    int occupTotal{1}, exact{1}, nuclExact{0}, occupElder{1}, step{1}, occupStep{}, ffo{};
    double isMut{0}, isOcc{0}, prOcc{0};

// Объявлеем вектор-матрицу
    vector<vector<int>> matrix;
    matrix.push_back(vector<int>());

// Объявлеем вектор-ДНК
    vector<int> occupNumer(nichNumb);
//* 3.1 -------------------------------------------




//* ANCHOR 3.2 -- ОСНОВНОЙ ЦИКЛИЗМ
    //^ 1 цикл - Повторения для получения доверительного интервала
    for (int conIntStep0 = conStepsZ; conIntStep0 <= conSteps2; conIntStep0++) {
        double jumpStep0 {jumpProbTy};

        //~ Создаём матрицу хранищую DNA для вида каждой ниши
        int **niches = matrxDyn(nichNumb, dnaSize);

        occupNumer.resize(nichNumb);
        occupNumer[0] = 0;
        std::copy(digiotide, digiotide + dnaSize, niches[0]);

        //~ Обнуляем все остальные ниши (заполняем их DNA четвёрками "4")
        for (int fours = 1; fours < nichNumb; fours++)
        {
            std::copy(digiotideFour, digiotideFour + dnaSize, niches[fours]);
            occupNumer[fours] = 0;
            /*
            std::cout << endl << fours << " ";
            for (int iii{0}; iii < dnaSize; iii++){
                std::cout << to_string(niches[fours][iii]);
            }
            std::cout << endl;
            */
        }

        //~ Заполняем матрицу нулями
        matrix.resize(nichNumb);

        for (int row = 0; row < nichNumb; row++){
            matrix[row].resize(nichNumb);
            for (int col = 0; col <= nichNumb; col++){
                matrix[row][col] = 0;
            }
        }

        //~ Сбрасываем все переменные перед основным циклом заполнения ниш
        occupTotal  = 1;    // Сколько занято ниш
        occupElder  = 1;    // Старшая заполенная ниша
        exact       = 1;    // Ниша из которой будет переход
        isMut       = 0;    // Случайное число вероятности мутации (во всей DNA)
        nuclExact   = 0;    // Случайный нуклеотид в котором произойдёт мутация, если мутация состоится
        step        = 1;    // Шаг (поколение)
        isOcc       = 0;    // Случайное число вероятности перехода
        prOcc       = 0;    // Совокупная вероятность перехода

        ffo = nichNumb - occupTotal; //^ Оставшиеся свободные ниши
        //int ffo2 {1000};             //^ Если нужны дополнительные шаги после заполнения всех ниш +223 +279 (строки поменять)

        //^ Основной цикл -> повторять пока есть свободные ниши
        while (ffo > 0)
        {
            //~ 1 ЧАСТЬ - Мутации
            for (int iii = 0; iii < occupTotal; iii++)
            {
                isMut = mTwist_real(0, 1);

                if (isMut <= mutProbTy)
                {
                    nuclExact = mTwist_int(0, dnaSize - 1);

                    do
                    {
                        d = mTwist_int(0, 3);
                    } while (niches[iii][nuclExact] == d);

                    niches[iii][nuclExact] = d;
                }
            }

            //~ 2 ЧАСТЬ - Заполнение свободных ниш
            occupStep  = occupTotal;        // Заполнено ниш на начало шага
            occupElder = occupTotal;        // Старшая заполненая ниша
            

            while (occupStep > 0 && ffo > 0)
            {
                isOcc = mTwist_real(0, 1);                    //^ Случайная вероятность (должна быть меньше prOcc)
                prOcc = 1 - pow((1 - jumpStep0), occupStep);  //^ Совокупная вероятность того что "хотя бы один вид" займёт новую нишу

                
                if (isOcc < prOcc) //^ Если условие выполняется, то проиходит переход в свободную нишу
                {
                    exact = mTwist_int(0, (occupElder - 1));    // Выбираем из какой ниши будет переход
                    std::copy(niches[exact], niches[exact] + dnaSize, niches[occupTotal]);

                    for (int iii = 0; iii < occupTotal; iii++)
                    {
                        if (iii < exact)
                            matrix[occupTotal][iii] = matrix[exact][iii];
                        else
                            matrix[occupTotal][iii] = matrix[iii][exact];
                    }

                    //std::cout << endl << "|| From " << exact+1 << " to " << occupTotal+1 << " || step " << step << endl;

                    matrix[occupTotal][exact] = step;
                    occupTotal++;
                    occupStep--;
                    ffo--;
                    continue;
                }
                else
                    break;
            }
            
            //if (ffo <= 0) ffo2--; else 
            step++;
        }
        //^ Все ниши заняты! Конец повторения.




//* ANCHOR 3.3 -- Создаём файлы для матриц и количества шагов
        string filepath = "results/" + to_str(conIntStep0) + "s_matrixResult.txt"; // путь к файлу с изменяемым именем
        ofstream outM(filepath);

        if (!outM)
        {
            std::cerr << endl
                      << to_string(conIntStep0) + "_matrixResult.txt couldn't be opened for writing!" << endl;
        }

        int steps{step - 1};

        string filepath2 = "results/All_steps.txt"; // путь к файлу с изменяемым именем
        ofstream outS(filepath2, std::ios::app);

        if (!outS)
        {
            std::cerr << endl
                        << "All_steps.txt couldn't be opened for writing!" << endl;
        }

        outS << conIntStep0 << " " << steps << endl;
        outS.close(); //Закрываем файл

        long matrix_dist[nichNumb][nichNumb];

        for (int row = 0; row < nichNumb; row++)
        {
            matrix_dist[row][row] = 0; // Для того чтобы диагональ матрицы содержала нули
            for (int col = 0; col < row; col++)
            {
                matrix_dist[row][col] = (steps - matrix[row][col]) * 2;
                outM << matrix_dist[row][col] << "\t";
            }
            outM << matrix_dist[row][row] << "\t";
            outM << endl;
        }

//* 3.3 --------------------------------------------




//* ANCHOR 3.4 -- Записываем последовательности DNA
        // Создаём файл
        string filepath1 = "results/" + to_string(conIntStep0) + "s_DNAResult.txt"; // путь к файлу с изменяемым именем
        ofstream outD(filepath1, std::ios::out);

        if (!outD)
        {
            std::cerr << endl
                      << to_string(conIntStep0) + "s_DNAResult.txt couldn't be opened for writing!" << endl;
        }
            
        // Записываем последовательности DNA
        for (int iii{0}; iii < nichNumb; iii++)
        {
            outD << ">sp_" << iii + 1 << endl;

            for (int jjj{0}; jjj < dnaSize; jjj++)
            {
                d = niches[iii][jjj];

                if (d < 2)
                {
                    if (d == 0)
                        outD << 'T';
                    else
                        outD << 'C';
                }
                else if (d == 2)
                    outD << 'A';
                else
                    outD << 'G';
            }
            outD << endl;
        }
        
        outD.close(); //Закрываем файл
//* 3.4 --------------------------------------------

        // Удаляем матрицу
        for (int count = 0; count < nichNumb; count++)
        {
            delete[] niches[count];
        }
        delete[] niches;
    
        std::cout << endl << "|| GLOBAL STEP ||" << conIntStep0 << " / " << conSteps2 << endl;
    }
//* 3.2 --------------------------------------------

    std::cout << endl
              << "|| All processes completed. ||" << endl;
//* 3 ----------------------------------------------

    // Переводим строку и задерживаем окно консоли
    std::cout << endl;
    stopClosing();

    return 0;
    }
//~------------------------------------------------------------------------------------------
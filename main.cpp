#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <bitset>
#include <codecvt>
#include <locale>

using namespace std;

// Функция для вычисления статистических характеристик
void calcStats(const vector<unsigned char>& data, int maxCombLen, int shiftACF) {
    int totalBits = data.size() * 8;
    cout << "Объем выборки: " << totalBits << " бит" << endl;

    // Вычисление вероятностей появления комбинаций длиной от 1 до 4 бит
    for (int combLen = 1; combLen <= maxCombLen; combLen++) {
        unordered_map<unsigned int, int> combFreqMap;
        int totalCombs = totalBits - combLen + 1;

        for (size_t i = 0; i < data.size(); i++) {
            for (int j = 0; j <= 8 - combLen; j++) {
                unsigned int comb = (data[i] >> (8 - combLen - j)) & ((1 << combLen) - 1);
                combFreqMap[comb]++;
            }
        }

        cout << "Вероятности появления комбинаций длиной " << combLen << " бит:" << endl;
        for (auto p : combFreqMap) {
            double prob = (double)p.second / totalCombs;
            cout << "Комбинация " << bitset<4>(p.first).to_string().substr(4 - combLen) << ": " << prob << endl;
        }
        cout << endl;
    }

    // Вычисление вероятностей появления комбинаций длиной от 2 до 4 бит с одинаковым количеством единиц
    for (int combLen = 2; combLen <= maxCombLen; combLen++) {
        unordered_map<unsigned int, int> combFreqMap;
        int totalCombs = totalBits - combLen + 1;

        for (size_t i = 0; i < data.size(); i++) {
            for (int j = 0; j <= 8 - combLen; j++) {
                unsigned int comb = (data[i] >> (8 - combLen - j)) & ((1 << combLen) - 1);
                int onesCount = bitset<4>(comb).count();
                if (onesCount == combLen / 2) {
                    combFreqMap[comb]++;
                }
            }
        }

        cout << "Вероятности появления комбинаций длиной " << combLen << " бит с " << combLen / 2 << " единицами:" << endl;
        for (auto p : combFreqMap) {
            double prob = (double)p.second / totalCombs;
            cout << "Комбинация " << bitset<4>(p.first).to_string().substr(4 - combLen) << ": " << prob << endl;
        }
        cout << endl;
    }

    // Вычисление АКФ
    vector<double> acf(shiftACF + 1, 0.0);
    double mean = 0;

    for (unsigned char byte : data) {
        for (int i = 0; i < 8; i++) {
            mean += (byte >> (7 - i)) & 0x1;
        }
    }
    mean /= totalBits;

    for (int shift = 0; shift <= shiftACF; shift++) {
        for (size_t i = 0; i < data.size(); i++) {
            for (int j = 0; j < 8; j++) {
                if (i * 8 + j + shift < static_cast<size_t>(totalBits)) {
                    int bit1 = (data[i] >> (7 - j)) & 0x1;
                    int bit2 = (data[(i * 8 + j + shift) / 8] >> (7 - ((i * 8 + j + shift) % 8))) & 0x1;
                    acf[shift] += (bit1 - mean) * (bit2 - mean);
                }
            }
        }
        acf[shift] /= (totalBits - shift);
    }

    cout << "Автокорреляционная функция ПСП:" << endl;
    for (int shift = 0; shift <= shiftACF; shift++) {
        cout << "Сдвиг " << shift << ": " << acf[shift] << endl;
    }
    cout << endl;
}

// Функция для генерации гаммы с преобладанием
void generateGammaWithBias(vector<unsigned char>& gamma, int& rlz, int polynom, int n) {
    int onesCount = 0;

    for (size_t i = 0; i < gamma.size(); i++) {
        gamma[i] = 0;
        for (int j = 0; j < 8; j++) {
            int bitGamma;

            if (polynom == 1) {
                bitGamma = ((rlz >> 4) ^ (rlz >> 1)) & 0x1;
                rlz = (rlz << 1) | bitGamma;
                rlz &= 0x1F;
            } else {
                bitGamma = ((rlz >> 5) ^ (rlz >> 4) ^ (rlz >> 1) ^ rlz) & 0x1;
                rlz = (rlz << 1) | bitGamma;
                rlz &= 0x3F;
            }

            if (bitGamma == 1) {
                onesCount++;
                if (onesCount % n == 1) {
                    bitGamma = 0;  // Инвертирование каждой n-й единицы, начиная со второй
                }
            }

            gamma[i] = (gamma[i] << 1) | bitGamma;
        }
    }
}

// Функция для перекодирования данных в UTF-16 LE
void convertToUTF16LE(const vector<unsigned char>& data, vector<char>& utf16LEData) {
    utf16LEData.resize(data.size() * 2);
    for (size_t i = 0; i < data.size(); i++) {
        unsigned short value = static_cast<unsigned short>(data[i]);
        utf16LEData[i * 2] = static_cast<char>(value & 0xFF);
        utf16LEData[i * 2 + 1] = static_cast<char>((value >> 8) & 0xFF);
    }
}

int main() {
    int rlz1 = 31;
    int rlz2 = 51;
    char inName[256];
    char outName[256];
    int maskMode;
    int biasN;
    int combLen;
    int shiftACF;
    int polynom;
    int useUTF16LE;

    cout << "========= Меню =========" << endl;
    cout << "1. Ввод имени входного файла" << endl;
    cout << "2. Ввод имени выходного файла" << endl;
    cout << "3. Выбор режима маскирования (0 - без преобладания, 1 - с преобладанием)" << endl;
    cout << "4. Ввод значения n для преобладания (если выбран режим с преобладанием)" << endl;
    cout << "5. Ввод длины комбинаций для вычисления вероятностей" << endl;
    cout << "6. Ввод значения сдвига для вычисления АКФ" << endl;
    cout << "7. Выбор полинома для генерации гаммы (1 - x^5+x^2+1, 2 - x^6+x^5+x^2+x+1)" << endl;
    cout << "8. Применить перекодировку в UTF-16 LE (0 - нет, 1 - да)" << endl;
    cout << "9. Запуск маскирования и анализа" << endl;

    cout << "Введите имя входного файла: ";
    cin >> inName;
    cout << "Введите имя выходного файла: ";
    cin >> outName;
    cout << "Выберите режим маскирования (0 - без преобладания, 1 - с преобладанием): ";
    cin >> maskMode;

    if (maskMode == 1) {
        cout << "Введите значение n для преобладания: ";
        cin >> biasN;
    }

    cout << "Введите длину комбинаций для вычисления вероятностей: ";
    cin >> combLen;
    cout << "Введите значение максимального сдвига tau для вычисления АКФ: ";
    cin >> shiftACF;
    cout << "Выберите полином для генерации гаммы (1 - x^5+x^2+1, 2 - x^6+x^5+x^2+x+1): ";
    cin >> polynom;
    cout << "Применить перекодировку в UTF-16 LE (0 - нет, 1 - да): ";
    cin >> useUTF16LE;

    ifstream inFile(inName, ios::binary);
    ofstream outFile(outName, ios::binary);

    vector<unsigned char> fileData((istreambuf_iterator<char>(inFile)), istreambuf_iterator<char>());

    cout << "Статистические характеристики исходного файла:" << endl;
    calcStats(fileData, combLen, shiftACF);

    vector<unsigned char> gamma(fileData.size() - (fileData.size() > 54 ? 54 : 0));

    if (maskMode == 0) {
        // Генерация гаммы без преобладания
        for (size_t i = 0; i < gamma.size(); i++) {
            gamma[i] = 0;
            for (int j = 0; j < 8; j++) {
                gamma[i] = gamma[i] << 1;
                int bitGamma;

                if (polynom == 1) {
                    bitGamma = ((rlz1 >> 4) ^ (rlz1 >> 1)) & 0x1;
                    rlz1 = (rlz1 << 1) | bitGamma;
                    rlz1 &= 0x1F;
                } else {
                    bitGamma = ((rlz2 >> 5) ^ (rlz2 >> 4) ^ (rlz2 >> 1) ^ rlz2) & 0x1;
                    rlz2 = (rlz2 << 1) | bitGamma;
                    rlz2 &= 0x3F;
                }

                gamma[i] = gamma[i] + bitGamma;
            }
        }
    } else {
        // Генерация гаммы с преобладанием
        if (polynom == 1) {
            generateGammaWithBias(gamma, rlz1, polynom, biasN);
        } else {
            generateGammaWithBias(gamma, rlz2, polynom, biasN);
        }
    }

    // Маскирование данных
    if (fileData.size() > 54) {
        // Сохранение заголовка BMP
        vector<unsigned char> bmpHeader(fileData.begin(), fileData.begin() + 54);

        // Маскирование файла BMP без заголовка
        for (size_t i = 54; i < fileData.size(); i++) {
            fileData[i] ^= gamma[i - 54];
        }

        // Запись сохраненного заголовка BMP в начало файла
        copy(bmpHeader.begin(), bmpHeader.end(), fileData.begin());
    } else {
        // Маскирование обычного файла
        for (size_t i = 0; i < fileData.size(); i++) {
            fileData[i] ^= gamma[i];
        }
    }

if (useUTF16LE) {
    // Перекодирование данных в UTF-16 LE
    vector<char> utf16LEData;
    convertToUTF16LE(fileData, utf16LEData);

    // Запись данных в выходной файл в кодировке UTF-16 LE
    outFile.write(utf16LEData.data(), utf16LEData.size());
} else {
    // Запись данных в выходной файл без перекодировки
    outFile.write(reinterpret_cast<const char*>(fileData.data()), fileData.size());
}

cout << "Статистические характеристики замаскированного файла:" << endl;
calcStats(fileData, combLen, shiftACF);

inFile.close();
outFile.close();

return 0;
}
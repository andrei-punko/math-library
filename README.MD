# Библиотека численных методов

[![C/C++ CI](https://github.com/andrei-punko/math-library/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/andrei-punko/math-library/actions/workflows/c-cpp.yml)

Решатели дифференциальных уравнений в частных производных (и сопутствующий код) мигрированы на Java и перенесены в другой [репозиторий](https://github.com/andrei-punko/pde-solvers).
Стоит смотреть его, тогда как текущий репозиторий представляет в основном археологическую ценность.

## Содержание

- Аппроксимация функции заданной таблично с заданным среднеквадратичным отклонением
  при помощи произвольного набора функций.

- Численные значения констант и параметров веществ

- Генерация равномерно распределенной в заданном интервале случайной величины.  
  Алгоритм взят из книги "Numerical Recipes in C", с.279

- Генерация нормально распределенной случайной величины с заданными мат. ожиданием и дисперсией.  
  Алгоритм взят из книги Кнут "Искусство программирования", раздел 3.4, с.4

- Для функции, заданной таблично реализованы:
  дифференцирование (1,2,3 производные), интегрирование, интерполяция, табулирование.
  Количество точек в задающей матрице может быть изменено.

- Метод прогонки

- Класс CInterval: интервал [x1; x2], шаг h на нем, число кусков N, на которые он разбивается шагом

- Класс CMatrix: одно- или двумерный массив (размерность = кол-ву переменных в конструкторе)

- Интегрирование методом Симпсона (интервал интегрирования может содержать нечетное количество частей)

- Табулирование функции на заданном интервале и сохранение результатов в файл. 3 варианта реализации

- Нахождение минимума (максимума) функции на заданном интервале методом золотого сечения

- Решение системы линейных уравнений методом Гаусса с выбором главного элемента

- Нахождение обратной матрицы

- Вычисление определителя матрицы (приведением ее к треугольной)

- Математическое округление

- Сохранение плотности распределения набора чисел с заданным шагом

- Сохранение плотности распределения функции на заданном интервале с заданным шагом

- Класс CLattice инкапсулирует параметры решетки частиц, расположенной в кубе:

      ребро куба и его половина            L, L_2(для ускорения вычислений)
      масса одной частицы                  m_1
      количество частиц в кубе             N
      парный потенциал и его производная   V, dV_dr(для ускорения вычислений)
      координаты и скорости частиц         XYZ, VXYZ

  После создания объекта класса вызовом одной из ф-й Init...() выбирается тип решетки:

      InitAmorph(N_v)   аморфное в-во; частицы с концентрацей N_v равномерно размещаются в кубе случайным образом
      InitVCC(a)        ОЦК решетка с параметром а
      InitManual(m)     координаты каждого из атомов задаются "вручную"

- Класс CUnSize реализует переход между знаковыми и беззнаковыми величинами. Масштаб задается набором единиц x,m,t

- Метод молекулярной динамики (не ручаюсь за его правильность)

- Решение ОДУ dU/dx = F(x,U) методом Рунге-Кутта 4 порядка

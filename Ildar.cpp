#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
typedef long long
    ll; // назначили короткое имя типу long long для удобства написания кода:
// вместо "long long x" можно теперь писать "ll x".
typedef long double ld; // аналогично ll
ld rand(ld minimal, ld maximal, bool include_minimal = true) {
  if (minimal > maximal)
    swap(minimal, maximal);
  ld r = rand();
  while (!include_minimal && r == 0)
    r = rand();
  r /= RAND_MAX;
  return minimal + (maximal - minimal) * r;
}
ld Umax = 100;
struct ModelVariables {
  ll preys, predators;
};
struct ModelVariablesf {
  ld preys, predators;
};
struct Results {
  vector<ld> t;
  vector<ModelVariables> u;
};
struct Resultsf {
  vector<ld> t;
  vector<ModelVariablesf> u;
};
Results gillespie_algorithm(ld Tmax, ld a, ld b, ld c, ld d, ll preys0,
                            ll predators0) {
  /*
  Реализуем стохастический вариант системы ОДУ вида:
  du/dt = a * u - b * u * v  -- динамика роста травы на поле: она растёт,
  её едят;
  dv/dt = c * u - d * v      -- динамика изменения популяции
  травоядных на поле: чем больше корма, тем больше особей приходит пастись, но
  часть со временем уходит;
  */
  Results r;
  // сохраним состояние системы на начальный момент времени
  r.t = {0};
  r.u = {{preys0, predators0}};
  ld t = 0;
  ll preys_now = preys0, predators_now = predators0;
  ld au, buv, cu, dv, total, t_new;
  // пока не достигнуто целевое время, будем совершать событие за событием
  while (t < Tmax) {
    ld r1 = rand(0, 1, false);
    ld r2 = rand(0, 1);
 
    // для стохастического алгоритма надо разобрать систему уравнений на
    // компоненты, описывающие скорости отдельных процессов
    au = a * preys_now*(1-preys_now/Umax); // скорость роста травы
    buv = b * preys_now * predators_now; // скорость поедания травы скотом
    cu = c * preys_now * predators_now;     // скорость прихода скота
    dv = d * predators_now; // скорость ухода скота
    // нам потребуется суммарная скорость процессов
    total = au + buv + cu + dv;
    if (total == 0) {
      // будем экстренно завершать расчёты, если все процессы остановились или
      // наоборот, слишком сильно разогнались
      break;
    }
    // время ближайшего события определеяется в некоторой степени случайно
    // согласно вероятностному распределению
    t_new = t + 1 / total * log(1 / r1);
    // на основании датчика случайных чисел определяем, какое из событий
    // произошло
    if (r2 < au / total) {
      // выросла травинка
      preys_now++;
    } else if (r2 < (au + buv) / total) {
      // травинку съели
      preys_now--;
    } else if (r2 < (au + buv + cu) / total) {
      // пришло парнокопытное
      predators_now++;
    } else {
      // ушло парнокопытное
      predators_now--;
    }
    // сохраняем новую временную метку и численности
    r.t.push_back(t = t_new);
    r.u.push_back({preys_now, predators_now});
  }
  return r;
}
Resultsf ode_algorithm(ld Tmax, ld a, ld b, ld c, ld d, ll preys0,
                       ll predators0, ld epsilon, ll NPointsSave) {
  /*
    Реализуем численное решение системы ОДУ вида:
    du/dt = a * u - b * u * v  -- динамика роста травы на поле: она растёт,
    её едят;
    dv/dt = c * u - d * v      -- динамика изменения популяции
    травоядных на поле: чем больше корма, тем больше особей приходит пастись, но
    часть со временем уходит;
  */
  auto solve = [&](ll N, ll N0) {
    ld dt = Tmax / N;
    ll NN0 = N / N0;
    Resultsf r;
    // сохраним состояние системы на начальный момент времени
    r.t = {0};
    r.u = {{ld(preys0), ld(predators0)}};
    
    ld d_preys, d_predators, preys_now = preys0, predators_now = predators0,
                             preys_pred, predators_pred;
    // лямбда-функция, которая посчитает нам производную для жертв
    auto dPreys = [&](ld preys, ld predators) {
      return a * preys*(1-preys/Umax) - b * preys * predators;
    };
    // лямбда-функция, которая посчитает нам производную для хищников
    auto dPredators = [&](ld preys, ld predators) {
      return c * preys * predators - d * predators;
    };
    for (ll i = 0; i < N; i++) {
      // используем для ускорения сходимости метод предиктор-корректор
      // вычисляем производную
      d_preys = dPreys(preys_now, predators_now);
      d_predators = dPredators(preys_now, predators_now);
      // делаем предварительный расчёт шага по времени
      preys_pred = preys_now + d_preys * dt;
      predators_pred = predators_now + d_predators * dt;
      // уточняем производную
      d_preys += dPreys(preys_pred, predators_pred);
      d_predators += dPredators(preys_pred, predators_pred);
      // вычисляем новую точку более точно
      preys_now = preys_now + d_preys * dt * 0.5;
      predators_now = predators_now + d_predators * dt * 0.5;
      // сохраняем точки с некоторой частотой для проверки точности и вывода
      // графиков
      if ((i + 1) % NN0 == 0) {
        r.t.push_back((i + 1) * dt);
        r.u.push_back({preys_now, predators_now});
      }
    }
    return r;
  };
  ll N = NPointsSave;
  auto solution_low = solve(N, NPointsSave),
       solution_high = solve(N * 2, NPointsSave);
  auto compare = [&]() {
    ld md = 0;
    ld diff_preys, diff_predators;
    ld low_preys, low_predators;
    ld high_preys, high_predators;
    for (ll i = 0; i < solution_low.t.size(); i++) {
      low_preys = solution_low.u[i].preys;
      high_preys = solution_high.u[i].preys;
      low_predators = solution_low.u[i].predators;
      high_predators = solution_high.u[i].predators;
      diff_preys = fabs(low_preys - high_preys);
      diff_predators = fabs(low_predators - high_predators);
      md = fmax(md, fmax(diff_preys, diff_predators));
    }
    return md;
  };
  ld cmp = compare();
  while (cmp > epsilon) {
    N *= 2;
    solution_low = solution_high;
    solution_high = solve(N * 2, NPointsSave);
    cmp = compare();
    printf("N = %lld, error = %0.16Lf.\n", N, cmp);
  }
  printf("Needed accuracy achieved with %lld steps.\n", N);
  return solution_high;
}
bool saveToPythonScript(Results r, Resultsf rf, string s) {
  ofstream fout(s);
  if (fout.is_open()) {
    fout << "#!/bin/python3"
         << endl; // header is needed to run PY file as Bash script
    // пишем в файл результаты расчётов алгоритма Гиллеспи
    fout << "t = [";
    for (ll i = 0; i < r.t.size(); i++)
      fout << (i ? "," : "") + to_string(r.t[i]);
    fout << "]" << endl << "preys = [";
    for (ll i = 0; i < r.u.size(); i++)
      fout << (i ? "," : "") + to_string(r.u[i].preys);
    fout << "]" << endl << "predators = [";
    for (ll i = 0; i < r.u.size(); i++)
      fout << (i ? "," : "") + to_string(r.u[i].predators);
    fout << "]" << endl;
    // пишем в файл результаты расчётов системы ОДУ
    fout << "ode_t = [";
    for (ll i = 0; i < rf.t.size(); i++)
      fout << (i ? "," : "") + to_string(rf.t[i]);
    fout << "]" << endl << "ode_preys = [";
    for (ll i = 0; i < rf.u.size(); i++)
      fout << (i ? "," : "") + to_string(rf.u[i].preys);
    fout << "]" << endl << "ode_predators = [";
    for (ll i = 0; i < rf.u.size(); i++)
      fout << (i ? "," : "") + to_string(rf.u[i].predators);
    fout << "]" << endl;
    // рисуем графики
    fout << "from matplotlib import pyplot as plt" << endl;
    fout << "plt.plot(t, preys, label = \"Grass\", color = \"g\")" << endl;
    fout << "plt.plot(t, predators, label = \"Cows\", color = \"r\")" << endl;
    fout << "plt.plot(ode_t, ode_preys, label = \"ODE Grass\", color = \"g\", "
            "ls = \"dashed\")"
         << endl;
    fout << "plt.plot(ode_t, ode_predators, label = \"ODE Cows\", color = "
            "\"r\", ls = \"dashed\")"
         << endl;
    fout << "plt.xlabel(\"Time\")" << endl;
    fout << "plt.ylabel(\"Count\")" << endl;
    fout << "plt.legend()" << endl;
    fout << "plt.grid()" << endl;
    fout << "plt.show()";
    fout.close();
    return true;
  }
  return false;
}
void runPythonScript(string s) { system(("python3 " + s).c_str()); }
int main() {
  srand(time(nullptr));
  ld a = 0.5;
  ld b = 0.01;
  ld c = 0.01;
  ld d = 0.15;
  
  string s = "example.py";
  auto r = gillespie_algorithm(100, a, b, c, d, 20, 1);
  auto rf = ode_algorithm(100, a, b, c, d, 20, 1, 1e-6, 100);
  if (saveToPythonScript(r, rf, s))
    runPythonScript(s);
  return 0;
}
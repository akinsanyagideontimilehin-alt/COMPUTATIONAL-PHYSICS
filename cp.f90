! ================================================================
! РЕШЕНИЕ УРАВНЕНИЯ ПЯТЬЮ МЕТОДАМИ
! SOLUTION OF EQUATION USING FIVE METHODS
! ================================================================

program solve_equation
    implicit none
    
    ! Декларации переменных / Variable declarations
    real :: a, b, c, x0, x1, x_next
    real :: root_dichotomy, root_newton, root_secant, root_iter, root_comb
    integer :: steps_dichotomy, steps_newton, steps_secant, steps_iter, steps_comb
    integer :: i, max_iter
    real :: eps, x_newton, x_iter
    real :: fa, fb, fc, fx, dfx, f_x0, f_x1
    real :: a_new, b_new
    real :: true_root
    
    ! Внешние функции / External functions
    real :: f, df, d2f, phi
    
    ! Параметры / Parameters
    eps = 1.0e-6
    max_iter = 1000
    
    write(*,*) '============================================================'
    write(*,*) 'РЕШЕНИЕ УРАВНЕНИЯ / SOLVING EQUATION:'
    write(*,*) 'f(x) = 0.04*(2x^3 - 5x^2 - 13x + 9) = 0 на [0, 2]'
    write(*,*) '============================================================'
    write(*,*)
    
    ! ============================================================
    ! 1. МЕТОД ДИХОТОМИИ (ДЕЛЕНИЯ ОТРЕЗКА ПОПОЛАМ)
    !    BISECTION METHOD (DICHOTOMY)
    ! ============================================================
    a = 0.0
    b = 2.0
    steps_dichotomy = 0
    
    if (f(a) * f(b) >= 0.0) then
        write(*,*) 'ПРЕДУПРЕЖДЕНИЕ: Нет смены знака на [0, 2]!'
        write(*,*) 'WARNING: No sign change in [0, 2]!'
        root_dichotomy = 0.0
    else
        do i = 1, max_iter
            steps_dichotomy = steps_dichotomy + 1
            c = (a + b) / 2.0
            
            ! Проверка сходимости / Convergence check
            if ((b - a) / 2.0 < eps .or. abs(f(c)) < eps) then
                root_dichotomy = c
                exit
            end if
            
            if (f(a) * f(c) < 0.0) then
                b = c
            else
                a = c
            end if
        end do
    end if
    
    write(*,*) '1. МЕТОД ДИХОТОМИИ / BISECTION METHOD:'
    write(*,*) '   Корень / Root = ', root_dichotomy, '  f(root) = ', f(root_dichotomy)
    write(*,*) '   Шагов / Steps = ', steps_dichotomy
    write(*,*)
    
    ! ============================================================
    ! 2. МЕТОД НЬЮТОНА (КАСАТЕЛЬНЫХ)
    !    NEWTON'S METHOD (TANGENT)
    ! ============================================================
    x_newton = 2.0  ! Начальное приближение / Initial guess
    steps_newton = 0
    root_newton = 0.0
    
    do i = 1, max_iter
        steps_newton = steps_newton + 1
        fx = f(x_newton)
        dfx = df(x_newton)
        
        if (abs(dfx) < 1.0e-12) then
            write(*,*) 'ОШИБКА: Производная слишком мала!'
            write(*,*) 'ERROR: Derivative too small!'
            exit
        end if
        
        x_next = x_newton - fx / dfx
        
        ! Проверка сходимости / Convergence check
        if (abs(x_next - x_newton) < eps .or. abs(f(x_next)) < eps) then
            root_newton = x_next
            exit
        end if
        
        x_newton = x_next
    end do
    
    write(*,*) '2. МЕТОД НЬЮТОНА / NEWTON''S METHOD:'
    write(*,*) '   Корень / Root = ', root_newton, '  f(root) = ', f(root_newton)
    write(*,*) '   Шагов / Steps = ', steps_newton
    write(*,*)
    
    ! ============================================================
    ! 3. МЕТОД ХОРД (СЕКУЩИХ)
    !    SECANT METHOD (CHORD)
    ! ============================================================
    x0 = 0.0
    x1 = 2.0
    steps_secant = 0
    root_secant = 0.0
    
    do i = 1, max_iter
        steps_secant = steps_secant + 1
        f_x0 = f(x0)
        f_x1 = f(x1)
        
        if (abs(f_x1 - f_x0) < 1.0e-12) then
            write(*,*) 'ОШИБКА: Деление на ноль!'
            write(*,*) 'ERROR: Division by zero!'
            exit
        end if
        
        x_next = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
        
        ! Проверка сходимости / Convergence check
        if (abs(x_next - x1) < eps .or. abs(f(x_next)) < eps) then
            root_secant = x_next
            exit
        end if
        
        x0 = x1
        x1 = x_next
    end do
    
    write(*,*) '3. МЕТОД ХОРД / SECANT METHOD (CHORD):'
    write(*,*) '   Корень / Root = ', root_secant, '  f(root) = ', f(root_secant)
    write(*,*) '   Шагов / Steps = ', steps_secant
    write(*,*)
    
    ! ============================================================
    ! 4. МЕТОД ИТЕРАЦИЙ (ПОСЛЕДОВАТЕЛЬНЫХ ПРИБЛИЖЕНИЙ)
    !    FIXED-POINT ITERATION METHOD
    ! ============================================================
    ! Приведем к виду x = phi(x): 13x = 2x^3 - 5x^2 + 9
    ! Transform to x = phi(x): 13x = 2x^3 - 5x^2 + 9
    ! => x = (2x^3 - 5x^2 + 9) / 13
    
    x_iter = 2.0  ! Начальное приближение / Initial guess
    steps_iter = 0
    root_iter = 0.0
    
    do i = 1, max_iter
        steps_iter = steps_iter + 1
        x_next = phi(x_iter)
        
        ! Проверка сходимости / Convergence check
        if (abs(x_next - x_iter) < eps .or. abs(f(x_next)) < eps) then
            root_iter = x_next
            exit
        end if
        
        x_iter = x_next
    end do
    
    write(*,*) '4. МЕТОД ИТЕРАЦИЙ / FIXED-POINT ITERATION:'
    write(*,*) '   Корень / Root = ', root_iter, '  f(root) = ', f(root_iter)
    write(*,*) '   Шагов / Steps = ', steps_iter
    write(*,*)
    
    ! ============================================================
    ! 5. КОМБИНИРОВАННЫЙ МЕТОД (ХОРД + КАСАТЕЛЬНЫХ)
    !    COMBINED METHOD (CHORD + NEWTON)
    ! ============================================================
    a = 0.0
    b = 2.0
    steps_comb = 0
    root_comb = 0.0
    
    if (f(a) * f(b) >= 0.0) then
        write(*,*) 'ПРЕДУПРЕЖДЕНИЕ: Нет смены знака на [0, 2]!'
        write(*,*) 'WARNING: No sign change in [0, 2]!'
    else
        do i = 1, max_iter
            steps_comb = steps_comb + 1
            
            ! Шаг хорд / Chord step
            if (abs(f(b) - f(a)) < 1.0e-12) then
                write(*,*) 'ОШИБКА: Деление на ноль!'
                write(*,*) 'ERROR: Division by zero!'
                exit
            end if
            a_new = a - f(a) * (b - a) / (f(b) - f(a))
            
            ! Шаг Ньютона / Newton step
            if (abs(df(b)) < 1.0e-12) then
                write(*,*) 'ОШИБКА: Производная слишком мала!'
                write(*,*) 'ERROR: Derivative too small!'
                exit
            end if
            b_new = b - f(b) / df(b)
            
            ! Обновление интервала / Update interval
            if (a_new < b_new) then
                a = a_new
                b = b_new
            else
                a = b_new
                b = a_new
            end if
            
            ! Проверка сходимости / Convergence check
            if (abs(b - a) < eps .or. (abs(f(a)) < eps .and. abs(f(b)) < eps)) then
                root_comb = (a + b) / 2.0
                exit
            end if
        end do
    end if
    
    write(*,*) '5. КОМБИНИРОВАННЫЙ МЕТОД / COMBINED METHOD (CHORD + NEWTON):'
    write(*,*) '   Корень / Root = ', root_comb, '  f(root) = ', f(root_comb)
    write(*,*) '   Шагов / Steps = ', steps_comb
    write(*,*)
    
    ! ============================================================
    ! СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ / SUMMARY TABLE
    ! ============================================================
    write(*,*) '============================================================'
    write(*,*) 'СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ / SUMMARY OF RESULTS:'
    write(*,*) '============================================================'
    write(*,*) 'Метод / Method                  Корень / Root     f(root)        Шагов / Steps'
    write(*,*) '----------------------------------------------------------------------'
    write(*,'(A25, F12.8, E12.2, I10)') '1. Дихотомия / Bisection', root_dichotomy, f(root_dichotomy), steps_dichotomy
    write(*,'(A25, F12.8, E12.2, I10)') '2. Ньютона / Newton', root_newton, f(root_newton), steps_newton
    write(*,'(A25, F12.8, E12.2, I10)') '3. Хорд / Secant', root_secant, f(root_secant), steps_secant
    write(*,'(A25, F12.8, E12.2, I10)') '4. Итераций / Fixed-point', root_iter, f(root_iter), steps_iter
    write(*,'(A25, F12.8, E12.2, I10)') '5. Комбинированный / Combined', root_comb, f(root_comb), steps_comb
    write(*,*) '============================================================'
    
    ! ============================================================
    ! ВЫЧИСЛЕНИЕ ТОЧНОГО КОРНЯ (ЭТАЛОН) / TRUE ROOT CALCULATION
    ! ============================================================
    ! Используем метод дихотомии с высокой точностью
    ! Using bisection with high precision
    a = 0.0
    b = 1.0
    do i = 1, 80
        c = (a + b) / 2.0
        if (f(c) > 0.0) then
            a = c
        else
            b = c
        end if
    end do
    true_root = (a + b) / 2.0
    
    write(*,*)
    write(*,*) 'ТОЧНЫЙ КОРЕНЬ / TRUE ROOT (80 итераций / 80 iterations):'
    write(*,*) '   x = ', true_root
    write(*,*) '   f(x) = ', f(true_root)
    write(*,*)
    
    write(*,*) '============================================================'
    write(*,*) 'ВСЕ МЕТОДЫ ДОЛЖНЫ СХОДИТЬСЯ К ОДНОМУ И ТОМУ ЖЕ КОРНЮ'
    write(*,*) 'ALL METHODS SHOULD CONVERGE TO THE SAME ROOT'
    write(*,*) '============================================================'
    
    ! ============================================================
    ! СОХРАНЕНИЕ РЕЗУЛЬТАТОВ В ФАЙЛ / SAVE RESULTS TO FILE
    ! ============================================================
    open(unit=10, file='results.txt', status='unknown')
    write(10,*) 'РЕЗУЛЬТАТЫ РЕШЕНИЯ УРАВНЕНИЯ / EQUATION SOLUTION RESULTS'
    write(10,*) '============================================================'
    write(10,*) 'f(x) = 0.04*(2x^3 - 5x^2 - 13x + 9) = 0  на [0, 2]'
    write(10,*) 'Точность / Precision: eps = ', eps
    write(10,*)
    write(10,*) 'Метод / Method                  Корень / Root     Шагов / Steps'
    write(10,*) '----------------------------------------------------------------'
    write(10,'(A25, F12.8, I10)') '1. Дихотомия / Bisection', root_dichotomy, steps_dichotomy
    write(10,'(A25, F12.8, I10)') '2. Ньютона / Newton', root_newton, steps_newton
    write(10,'(A25, F12.8, I10)') '3. Хорд / Secant', root_secant, steps_secant
    write(10,'(A25, F12.8, I10)') '4. Итераций / Fixed-point', root_iter, steps_iter
    write(10,'(A25, F12.8, I10)') '5. Комбинированный / Combined', root_comb, steps_comb
    write(10,*)
    write(10,*) 'Точный корень / True root: ', true_root
    close(10)
    
    write(*,*) 'РЕЗУЛЬТАТЫ СОХРАНЕНЫ В ФАЙЛ / RESULTS SAVED TO FILE: results.txt'
    
end program solve_equation

! ================================================================
! ФУНКЦИИ / FUNCTIONS
! ================================================================

! Исходная функция / Original function
function f(x)
    implicit none
    real :: f
    real, intent(in) :: x
    
    f = 0.04 * (2.0 * x**3 - 5.0 * x**2 - 13.0 * x + 9.0)
end function f

! Первая производная / First derivative
function df(x)
    implicit none
    real :: df
    real, intent(in) :: x
    
    df = 0.04 * (6.0 * x**2 - 10.0 * x - 13.0)
end function df

! Вторая производная / Second derivative
function d2f(x)
    implicit none
    real :: d2f
    real, intent(in) :: x
    
    d2f = 0.04 * (12.0 * x - 10.0)
end function d2f

! Функция для метода итераций / Function for iteration method
! x = phi(x) = (2x^3 - 5x^2 + 9) / 13
function phi(x)
    implicit none
    real :: phi
    real, intent(in) :: x
    
    phi = (2.0 * x**3 - 5.0 * x**2 + 9.0) / 13.0
end function phi
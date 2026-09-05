program root_finding
    implicit none

    ! Global constraints [Глобальные ограничения]
    double precision, parameter :: eps = 1d-6
    integer, parameter :: max_iter = 1000

    ! Bisection variables [Переменные для метода дихотомии]
    double precision :: a_bis, b_bis, c_bis, root_dichotomy
    integer :: steps_dichotomy

    ! Newton variables [Переменные для метода Ньютона]
    double precision :: x_newton, x_next_n, root_newton
    integer :: steps_newton, i

    ! Secant variables [Переменные для метода хорд]
    double precision :: x0_sec, x1_sec, x_next_s, root_secant, f_x0, f_x1
    integer :: steps_secant

    ! Fixed-point variables [Переменные для метода простых итераций]
    double precision :: x_iter, x_next_i, root_iter
    integer :: steps_iter

    ! Combined variables [Переменные для комбинированного метода]
    double precision :: a_comb, b_comb, a_new, b_new, root_comb
    integer :: steps_comb

    print *, "======================================================================"
    print *, "SOLVING: f(x) = 0.04*(2x^3 - 5x^2 - 13x + 9) = 0 on"
    print *, "======================================================================"
    print *, ""

    ! ============================================================
    ! 1. BISECTION METHOD (Dichotomy) [1. МЕТОД ДИХОТОМИИ (деления пополам)]
    ! ============================================================
    a_bis = 0.0d0
    b_bis = 2.0d0
    steps_dichotomy = 0
    root_dichotomy = 0.0d0

    if (f(a_bis) * f(b_bis) >= 0.0d0) then
        print *, "WARNING: No sign change in [0, 2] for bisection!"
    else
        do while (((b_bis - a_bis) / 2.0d0 > eps) .and. (steps_dichotomy < max_iter))
            steps_dichotomy = steps_dichotomy + 1
            c_bis = (a_bis + b_bis) / 2.0d0
            if (abs(f(c_bis)) < eps) then
                a_bis = c_bis
                b_bis = c_bis
                exit
            else if (f(a_bis) * f(c_bis) < 0.0d0) then
                b_bis = c_bis
            else
                a_bis = c_bis
            end if
        end do
        root_dichotomy = (a_bis + b_bis) / 2.0d0
    end if

    print '(A, F10.8, A, E10.2, A, I4)', &
        "1. BISECTION:    Root = ", root_dichotomy, &
        ", f(root) = ", f(root_dichotomy), ", Steps = ", steps_dichotomy

    ! ============================================================
    ! 2. NEWTON'S METHOD (Tangent) [2. МЕТОД НЬЮТОНА (касательных)]
    ! ============================================================
    x_newton = 2.0d0  ! Initial guess [Начальное приближение]
    steps_newton = 0
    root_newton = 0.0d0

    do i = 1, max_iter
        steps_newton = steps_newton + 1
        if (abs(df(x_newton)) < 1d-12) then
            print *, "   ERROR: Derivative too small in Newton method"
            exit
        end if
        x_next_n = x_newton - f(x_newton) / df(x_newton)
        if ((abs(x_next_n - x_newton) < eps) .or. (abs(f(x_next_n)) < eps)) then
            root_newton = x_next_n
            exit
        end if
        x_newton = x_next_n
    end do

    print '(A, F10.8, A, E10.2, A, I4)', &
        "2. NEWTON:       Root = ", root_newton, &
        ", f(root) = ", f(root_newton), ", Steps = ", steps_newton

    ! ============================================================
    ! 3. SECANT METHOD (Chord) [3. МЕТОД ХОРД (секущих)]
    ! ============================================================
    x0_sec = 0.0d0
    x1_sec = 2.0d0
    steps_secant = 0
    root_secant = 0.0d0

    do i = 1, max_iter
        steps_secant = steps_secant + 1
        f_x0 = f(x0_sec)
        f_x1 = f(x1_sec)
        if (abs(f_x1 - f_x0) < 1d-12) then
            print *, "   ERROR: Division by zero in secant method"
            exit
        end if
        x_next_s = x1_sec - f_x1 * (x1_sec - x0_sec) / (f_x1 - f_x0)
        if ((abs(x_next_s - x1_sec) < eps) .or. (abs(f(x_next_s)) < eps)) then
            root_secant = x_next_s
            exit
        end if
        x0_sec = x1_sec
        x1_sec = x_next_s
    end do

    print '(A, F10.8, A, E10.2, A, I4)', &
        "3. SECANT:       Root = ", root_secant, &
        ", f(root) = ", f(root_secant), ", Steps = ", steps_secant

    ! ============================================================
    ! 4. FIXED-POINT ITERATION [4. МЕТОД ПРОСТЫХ ИТЕРАЦИЙ]
    ! ============================================================
    x_iter = 2.0d0  ! Starting point [Начальная точка]
    steps_iter = 0
    root_iter = 0.0d0

    do i = 1, max_iter
        steps_iter = steps_iter + 1
        x_next_i = phi(x_iter)
        if ((abs(x_next_i - x_iter) < eps) .or. (abs(f(x_next_i)) < eps)) then
            root_iter = x_next_i
            exit
        end if
        x_iter = x_next_i
    end do

    print '(A, F10.8, A, E10.2, A, I4)', &
        "4. FIXED-POINT:  Root = ", root_iter, &
        ", f(root) = ", f(root_iter), ", Steps = ", steps_iter

    ! ============================================================
    ! 5. COMBINED METHOD [5. КОМБИНИРОВАННЫЙ МЕТОД]
    ! ============================================================
    a_comb = 0.0d0
    b_comb = 2.0d0
    steps_comb = 0
    root_comb = 0.0d0

    if (f(a_comb) * f(b_comb) >= 0.0d0) then
        print *, "WARNING: No sign change for combined method!"
    else
        do i = 1, max_iter
            steps_comb = steps_comb + 1
            
            ! Secant step from 'a' side [Шаг метода хорд со стороны 'a']
            if (abs(f(b_comb) - f(a_comb)) < 1d-12) exit
            a_new = a_comb - f(a_comb) * (b_comb - a_comb) / (f(b_comb) - f(a_comb))
            
            ! Newton step from 'b' side [Шаг метода Ньютона со стороны 'b']
            if (abs(df(b_comb)) < 1d-12) exit
            b_new = b_comb - f(b_comb) / df(b_comb)
            
            ! Maintain the boundaries [Сохранение границ интервала]
            if (a_new < b_new) then
                a_comb = a_new
                b_comb = b_new
            else
                a_comb = b_new
                b_comb = a_new
            end if
            
            ! Check convergence conditions [Проверка условий сходимости]
            if ((abs(b_comb - a_comb) < eps) .or. &
                ((abs(f(a_comb)) < eps) .and. (abs(f(b_comb)) < eps))) then
                root_comb = (a_comb + b_comb) / 2.0d0
                exit
            end if
        end do
    end if

    print '(A, F10.8, A, E10.2, A, I4)', &
        "5. COMBINED:     Root = ", root_comb, &
        ", f(root) = ", f(root_comb), ", Steps = ", steps_comb
    print *, ""

contains

    ! The main cubic function f(x) [Основная кубическая функция f(x)]
    double precision function f(x)
        double precision, intent(in) :: x
        f = 0.04d0 * (2.0d0 * x**3 - 5.0d0 * x**2 - 13.0d0 * x + 9.0d0)
    end function f

    ! The first derivative f'(x) [Первая производная f'(x)]
    double precision function df(x)
        double precision, intent(in) :: x
        df = 0.04d0 * (6.0d0 * x**2 - 10.0d0 * x - 13.0d0)
    end function df

    ! Fixed-point iteration function phi(x) [Итерационная функция phi(x)]
    double precision function phi(x)
        double precision, intent(in) :: x
        phi = (2.0d0 * x**3 - 5.0d0 * x**2 + 9.0d0) / 13.0d0
    end function phi

end program root_finding

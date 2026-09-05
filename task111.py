import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font for better readability
rcParams['font.size'] = 10
rcParams['axes.titlesize'] = 12

# ============================================================
# FUNCTION DEFINITIONS
# ============================================================
def f(x):
    return 0.04 * (2 * x**3 - 5 * x**2 - 13 * x + 9)

def df(x):
    return 0.04 * (6 * x**2 - 10 * x - 13)

def d2f(x):
    return 0.04 * (12 * x - 10)

def phi(x):
    """Fixed-point iteration function: x = phi(x)"""
    return (2 * x**3 - 5 * x**2 + 9) / 13

eps = 1e-6
max_iter = 1000

print("=" * 70)
print("SOLVING: f(x) = 0.04*(2x^3 - 5x^2 - 13x + 9) = 0 on [0, 2]")
print("=" * 70)
print()

# ============================================================
# 1. BISECTION METHOD (Dichotomy)
# ============================================================
a, b = 0.0, 2.0
bisection_history = []
if f(a) * f(b) >= 0:
    print("WARNING: No sign change in [0, 2] for bisection!")
    root_dichotomy = None
    steps_dichotomy = 0
else:
    steps_dichotomy = 0
    while (b - a) / 2 > eps and steps_dichotomy < max_iter:
        steps_dichotomy += 1
        c = (a + b) / 2
        bisection_history.append(c)
        if abs(f(c)) < eps:
            a = b = c
            break
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    root_dichotomy = (a + b) / 2

print(f"1. BISECTION METHOD:")
print(f"   Root = {root_dichotomy:.8f}, f(root) = {f(root_dichotomy):.2e}, Steps = {steps_dichotomy}")
print()

# ============================================================
# 2. NEWTON'S METHOD (Tangent)
# ============================================================
x_newton = 2.0
newton_history = [x_newton]
steps_newton = 0
root_newton = None
for i in range(max_iter):
    steps_newton += 1
    fx = f(x_newton)
    dfx = df(x_newton)
    if abs(dfx) < 1e-12:
        print("   ERROR: Derivative too small")
        break
    x_next = x_newton - fx / dfx
    newton_history.append(x_next)
    if abs(x_next - x_newton) < eps or abs(f(x_next)) < eps:
        root_newton = x_next
        break
    x_newton = x_next

print(f"2. NEWTON'S METHOD:")
print(f"   Root = {root_newton:.8f}, f(root) = {f(root_newton):.2e}, Steps = {steps_newton}")
print()

# ============================================================
# 3. SECANT METHOD (Chord)
# ============================================================
x0, x1 = 0.0, 2.0
secant_history = [x0, x1]
steps_secant = 0
root_secant = None
for i in range(max_iter):
    steps_secant += 1
    f_x0, f_x1 = f(x0), f(x1)
    if abs(f_x1 - f_x0) < 1e-12:
        print("   ERROR: Division by zero in secant method")
        break
    x_next = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
    secant_history.append(x_next)
    if abs(x_next - x1) < eps or abs(f(x_next)) < eps:
        root_secant = x_next
        break
    x0, x1 = x1, x_next

print(f"3. SECANT METHOD (Chord):")
print(f"   Root = {root_secant:.8f}, f(root) = {f(root_secant):.2e}, Steps = {steps_secant}")
print()

# ============================================================
# 4. FIXED-POINT ITERATION
# ============================================================
x_iter = 2.0
iteration_history = [x_iter]
steps_iter = 0
root_iter = None
for i in range(max_iter):
    steps_iter += 1
    x_next = phi(x_iter)
    iteration_history.append(x_next)
    if abs(x_next - x_iter) < eps or abs(f(x_next)) < eps:
        root_iter = x_next
        break
    x_iter = x_next

print(f"4. FIXED-POINT ITERATION:")
print(f"   Root = {root_iter:.8f}, f(root) = {f(root_iter):.2e}, Steps = {steps_iter}")
print()

# ============================================================
# 5. COMBINED METHOD (Chord + Newton)
# ============================================================
a, b = 0.0, 2.0
combined_history = [(a, b)]
steps_comb = 0
root_comb = None

if f(a) * f(b) >= 0:
    print("WARNING: No sign change in [0, 2] for combined method!")
else:
    for i in range(max_iter):
        steps_comb += 1
        
        # Secant step from a side
        if abs(f(b) - f(a)) < 1e-12:
            print("   ERROR: Division by zero in combined method")
            break
        a_new = a - f(a) * (b - a) / (f(b) - f(a))
        
        # Newton step from b side
        if abs(df(b)) < 1e-12:
            print("   ERROR: Derivative too small in combined method")
            break
        b_new = b - f(b) / df(b)
        
        # Update: keep the bracket
        if a_new < b_new:
            a, b = a_new, b_new
        else:
            a, b = b_new, a_new
        
        combined_history.append((a, b))
        
        if abs(b - a) < eps or (abs(f(a)) < eps and abs(f(b)) < eps):
            root_comb = (a + b) / 2
            break

print(f"5. COMBINED METHOD (Chord + Newton):")
print(f"   Root = {root_comb:.8f}, f(root) = {f(root_comb):.2e}, Steps = {steps_comb}")
print()

# ============================================================
# SUMMARY
# ============================================================
print("=" * 70)
print("SUMMARY OF RESULTS:")
print("=" * 70)
print(f"{'Method':<25} {'Root':<15} {'f(root)':<15} {'Steps':<10}")
print("-" * 70)
print(f"{'1. Bisection':<25} {root_dichotomy:<15.8f} {f(root_dichotomy):<15.2e} {steps_dichotomy:<10}")
print(f"{'2. Newton':<25} {root_newton:<15.8f} {f(root_newton):<15.2e} {steps_newton:<10}")
print(f"{'3. Secant':<25} {root_secant:<15.8f} {f(root_secant):<15.2e} {steps_secant:<10}")
print(f"{'4. Fixed-point':<25} {root_iter:<15.8f} {f(root_iter):<15.2e} {steps_iter:<10}")
print(f"{'5. Combined':<25} {root_comb:<15.8f} {f(root_comb):<15.2e} {steps_comb:<10}")
print()

# ============================================================
# GRAPH 1: FUNCTION PLOT WITH ROOT
# ============================================================
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

# --- Plot 1: Function f(x) ---
x_plot = np.linspace(-0.2, 2.2, 1000)
y_plot = f(x_plot)

ax1.plot(x_plot, y_plot, 'b-', linewidth=2, label='f(x)')
ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.8)

# Highlight the root
ax1.scatter(root_dichotomy, f(root_dichotomy), color='red', s=100, zorder=5, label=f'Root ≈ {root_dichotomy:.6f}')
ax1.axvline(x=root_dichotomy, color='red', linestyle='--', alpha=0.5)

# Highlight the interval
ax1.axvline(x=0, color='green', linestyle='--', alpha=0.5, label='Interval [0, 2]')
ax1.axvline(x=2, color='green', linestyle='--', alpha=0.5)

ax1.fill_between([0, 2], -0.5, 0.5, alpha=0.1, color='green')
ax1.set_xlabel('x')
ax1.set_ylabel('f(x)')
ax1.set_title('Function f(x) with Root')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_xlim(-0.2, 2.2)
ax1.set_ylim(-1.0, 0.5)

# --- Plot 2: Function zoomed near root ---
x_zoom = np.linspace(root_dichotomy - 0.1, root_dichotomy + 0.1, 500)
y_zoom = f(x_zoom)

ax2.plot(x_zoom, y_zoom, 'b-', linewidth=2, label='f(x)')
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax2.axvline(x=root_dichotomy, color='red', linestyle='--', alpha=0.7, label=f'Root = {root_dichotomy:.8f}')
ax2.scatter(root_dichotomy, f(root_dichotomy), color='red', s=100, zorder=5)

# Add a second root check
ax2.fill_between(x_zoom, -0.01, 0.01, alpha=0.1, color='yellow')
ax2.set_xlabel('x')
ax2.set_ylabel('f(x)')
ax2.set_title('Zoom near Root (x ∈ [r-0.1, r+0.1])')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_ylim(-0.02, 0.02)

# --- Plot 3: Convergence of methods ---
# Create x-axis for iterations
methods = ['Bisection', 'Newton', 'Secant', 'Fixed-point', 'Combined']
steps = [steps_dichotomy, steps_newton, steps_secant, steps_iter, steps_comb]
roots = [root_dichotomy, root_newton, root_secant, root_iter, root_comb]

# Calculate errors relative to the root
true_root = root_dichotomy  # Use bisection result as reference
errors = [abs(r - true_root) for r in roots]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

bars = ax3.bar(methods, steps, color=colors, alpha=0.7, edgecolor='black')

# Add step count on top of bars
for bar, step in zip(bars, steps):
    height = bar.get_height()
    ax3.text(bar.get_x() + bar.get_width()/2., height + 0.1,
             f'{step}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax3.set_xlabel('Method')
ax3.set_ylabel('Number of Iterations')
ax3.set_title('Iteration Count Comparison')
ax3.grid(True, alpha=0.3, axis='y')
ax3.set_ylim(0, max(steps) * 1.2)

plt.tight_layout()
plt.savefig('function_and_methods.png', dpi=300, bbox_inches='tight')
plt.show()

# ============================================================
# GRAPH 2: CONVERGENCE PLOTS FOR EACH METHOD
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

# 1. Bisection convergence
if bisection_history:
    ax = axes[0]
    ax.plot(range(1, len(bisection_history) + 1), bisection_history, 'b-o', markersize=4, linewidth=2)
    ax.axhline(y=true_root, color='red', linestyle='--', label=f'True root = {true_root:.6f}')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('x')
    ax.set_title('Bisection Method Convergence')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_xlim(0, len(bisection_history) + 1)

# 2. Newton convergence
if newton_history:
    ax = axes[1]
    ax.plot(range(len(newton_history)), newton_history, 'r-o', markersize=4, linewidth=2)
    ax.axhline(y=true_root, color='red', linestyle='--', label=f'True root = {true_root:.6f}')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('x')
    ax.set_title('Newton\'s Method Convergence')
    ax.grid(True, alpha=0.3)
    ax.legend()

# 3. Secant convergence
if secant_history:
    ax = axes[2]
    ax.plot(range(len(secant_history)), secant_history, 'g-o', markersize=4, linewidth=2)
    ax.axhline(y=true_root, color='red', linestyle='--', label=f'True root = {true_root:.6f}')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('x')
    ax.set_title('Secant Method Convergence')
    ax.grid(True, alpha=0.3)
    ax.legend()

# 4. Fixed-point convergence
if iteration_history:
    ax = axes[3]
    ax.plot(range(len(iteration_history)), iteration_history, 'm-o', markersize=4, linewidth=2)
    ax.axhline(y=true_root, color='red', linestyle='--', label=f'True root = {true_root:.6f}')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('x')
    ax.set_title('Fixed-Point Iteration Convergence')
    ax.grid(True, alpha=0.3)
    ax.legend()

# 5. Combined method convergence
if combined_history:
    ax = axes[4]
    a_hist = [h[0] for h in combined_history]
    b_hist = [h[1] for h in combined_history]
    ax.plot(range(len(a_hist)), a_hist, 'c-o', markersize=4, linewidth=2, label='a (secant)')
    ax.plot(range(len(b_hist)), b_hist, 'orange', marker='s', markersize=4, linewidth=2, label='b (Newton)')
    ax.axhline(y=true_root, color='red', linestyle='--', label=f'True root = {true_root:.6f}')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('x')
    ax.set_title('Combined Method Convergence')
    ax.grid(True, alpha=0.3)
    ax.legend()

# 6. Error comparison
ax = axes[5]
for method, hist, color, label in [
    ('Bisection', bisection_history, '#1f77b4', 'Bisection'),
    ('Newton', newton_history, '#ff7f0e', 'Newton'),
    ('Secant', secant_history, '#2ca02c', 'Secant'),
    ('Fixed-point', iteration_history, '#d62728', 'Fixed-point')
]:
    if hist:
        errors_hist = [abs(x - true_root) for x in hist]
        ax.semilogy(range(len(errors_hist)), errors_hist, 'o-', color=color, label=label, markersize=3, linewidth=1.5)

ax.set_xlabel('Iteration')
ax.set_ylabel('Error (log scale)')
ax.set_title('Error Comparison (log scale)')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_ylim(1e-12, 1e-1)

plt.tight_layout()
plt.savefig('convergence_plots.png', dpi=300, bbox_inches='tight')
plt.show()

# ============================================================
# GRAPH 3: VISUALIZATION OF EACH METHOD
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

x_plot = np.linspace(-0.5, 2.5, 500)
y_plot = f(x_plot)

# 1. Bisection visualization
ax = axes[0]
ax.plot(x_plot, y_plot, 'b-', linewidth=2)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
if bisection_history:
    # Show the first few bisection steps
    for i, c in enumerate(bisection_history[:10]):
        ax.axvline(x=c, color='gray', linestyle=':', alpha=0.5, linewidth=1)
        ax.scatter(c, f(c), color='red' if i == 0 else 'orange', s=30, alpha=0.7)
ax.scatter(true_root, f(true_root), color='red', s=100, zorder=5, label=f'Root = {true_root:.6f}')
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_title('Bisection Method: Interval Halving')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim(-0.2, 2.2)

# 2. Newton visualization
ax = axes[1]
ax.plot(x_plot, y_plot, 'b-', linewidth=2)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
if newton_history:
    # Show tangent lines for first few iterations
    for i in range(min(3, len(newton_history) - 1)):
        x0 = newton_history[i]
        y0 = f(x0)
        slope = df(x0)
        x_tangent = np.linspace(x0 - 0.5, x0 + 0.5, 100)
        y_tangent = y0 + slope * (x_tangent - x0)
        ax.plot(x_tangent, y_tangent, '--', alpha=0.5, linewidth=1)
        ax.scatter(x0, y0, color='orange', s=50)
ax.scatter(true_root, f(true_root), color='red', s=100, zorder=5, label=f'Root = {true_root:.6f}')
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_title('Newton\'s Method: Tangent Lines')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim(-0.2, 2.2)

# 3. Secant visualization
ax = axes[2]
ax.plot(x_plot, y_plot, 'b-', linewidth=2)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
if secant_history:
    # Show secant lines for first few iterations
    for i in range(min(3, len(secant_history) - 2)):
        x0, x1 = secant_history[i], secant_history[i+1]
        y0, y1 = f(x0), f(x1)
        if abs(y1 - y0) > 1e-12:
            x_secant = np.linspace(min(x0, x1) - 0.3, max(x0, x1) + 0.3, 100)
            slope = (y1 - y0) / (x1 - x0)
            y_secant = y0 + slope * (x_secant - x0)
            ax.plot(x_secant, y_secant, '--', alpha=0.5, linewidth=1)
        ax.scatter(x0, y0, color='orange', s=50)
ax.scatter(true_root, f(true_root), color='red', s=100, zorder=5, label=f'Root = {true_root:.6f}')
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_title('Secant Method: Secant Lines')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim(-0.2, 2.2)

# 4. Fixed-point visualization
ax = axes[3]
# Plot y = phi(x) and y = x
x_fixed = np.linspace(-0.2, 2.2, 500)
y_phi = phi(x_fixed)
ax.plot(x_fixed, y_phi, 'g-', linewidth=2, label='y = φ(x)')
ax.plot(x_fixed, x_fixed, 'r--', linewidth=1.5, label='y = x')
if iteration_history:
    # Show the iterative steps
    for i in range(min(5, len(iteration_history) - 1)):
        x0 = iteration_history[i]
        x1 = iteration_history[i+1]
        # Vertical line from (x0, x0) to (x0, φ(x0))
        ax.plot([x0, x0], [x0, phi(x0)], 'gray', linestyle=':', alpha=0.5)
        # Horizontal line from (x0, φ(x0)) to (φ(x0), φ(x0))
        ax.plot([x0, x1], [phi(x0), phi(x0)], 'gray', linestyle=':', alpha=0.5)
        ax.scatter(x0, phi(x0), color='orange', s=50)
ax.scatter(true_root, true_root, color='red', s=100, zorder=5, label=f'Root = {true_root:.6f}')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Fixed-Point Iteration: x = φ(x)')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim(-0.2, 2.2)
ax.set_ylim(-0.2, 2.2)

# 5. Combined method visualization
ax = axes[4]
ax.plot(x_plot, y_plot, 'b-', linewidth=2)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
if combined_history:
    # Show the brackets
    for i, (a, b) in enumerate(combined_history[:8]):
        ax.plot([a, b], [f(a), f(b)], 'orange', alpha=0.3, linewidth=1)
        ax.scatter(a, f(a), color='orange', s=20, alpha=0.5)
        ax.scatter(b, f(b), color='orange', s=20, alpha=0.5)
ax.scatter(true_root, f(true_root), color='red', s=100, zorder=5, label=f'Root = {true_root:.6f}')
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_title('Combined Method: Chord + Newton')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_xlim(-0.2, 2.2)

# 6. Method comparison overview
ax = axes[5]
# Plot the function with all method results
ax.plot(x_plot, y_plot, 'b-', linewidth=2)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)

# Mark results from each method
markers = ['o', 's', '^', 'D', '*']
method_names = ['Bisection', 'Newton', 'Secant', 'Fixed-point', 'Combined']
method_roots = [root_dichotomy, root_newton, root_secant, root_iter, root_comb]
method_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

for i, (name, root, color, marker) in enumerate(zip(method_names, method_roots, method_colors, markers)):
    ax.scatter(root, f(root), color=color, marker=marker, s=150, zorder=5, 
               label=f'{name}: {root:.6f}', edgecolor='black', linewidth=1)

ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_title('All Methods: Roots Comparison')
ax.grid(True, alpha=0.3)
ax.legend(loc='lower left', fontsize=8)
ax.set_xlim(-0.2, 2.2)
ax.set_ylim(-0.5, 0.5)

plt.tight_layout()
plt.savefig('methods_visualization.png', dpi=300, bbox_inches='tight')
plt.show()

print("=" * 70)
print("GRAPHS SAVED:")
print("  1. function_and_methods.png - Function plot and iteration comparison")
print("  2. convergence_plots.png - Convergence of each method")
print("  3. methods_visualization.png - Visual explanation of each method")
print("=" * 70)
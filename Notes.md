
- [1. Code](#1-code)
  - [1.1. run\_dense3D\_magfield](#11-run_dense3d_magfield)
    - [1.1.1. pipeline](#111-pipeline)
    - [1.1.2. dynModel](#112-dynmodel)
    - [1.1.3. domain\_cartesian\_dx](#113-domain_cartesian_dx)
    - [1.1.4. particleFilter](#114-particlefilter)
      - [1.1.4.1. ai](#1141-ai)
      - [1.1.4.2. xl](#1142-xl)
      - [1.1.4.3. particle](#1143-particle)
      - [1.1.4.4. how to use GP to predict](#1144-how-to-use-gp-to-predict)
      - [1.1.4.5. sequential process](#1145-sequential-process)
- [2. Q\&A](#2-qa)
  - [2.1. why does dynModel contain noise when predicting state?](#21-why-does-dynmodel-contain-noise-when-predicting-state)
  - [2.2. SS = dy\* P(:,:,i)\*dy' + R;](#22-ss--dy-pidy--r)
- [3. Issue](#3-issue)
- [4. Todo](#4-todo)
  - [4.1. compare smoother with filter](#41-compare-smoother-with-filter)
  - [4.2. try to understand the dynModel: does it contain anything related to angular velocity?](#42-try-to-understand-the-dynmodel-does-it-contain-anything-related-to-angular-velocity)
  - [4.3. 3.3](#43-33)

# 1. Code
## 1.1. run_dense3D_magfield
### 1.1.1. pipeline
![alt text](figs/ancestor.png)
### 1.1.2. dynModel
```
    dynModel = @(xn,dx,dt,Q) [xn(1:2) + [cos(xn(3)), -sin(xn(3)) ; ...
                sin(xn(3)), cos(xn(3))]' * dx(1:2)' ; xn(3) + dx(3) + sqrt(Q) * randn]; 
```
### 1.1.3. domain_cartesian_dx

```
[eigenval,~,eigenfun_dx,NN] = domain_cartesian_dx(nBasisFunctions,d,LL);
```

```
domain_cartesian_dx - Laplace operator eigendecomposition in a hypercube
```
![alt text](figs/bxfun.png)

```
  eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
```


### 1.1.4. particleFilter
#### 1.1.4.1. ai
a means ancester
#### 1.1.4.2. xl
xl = xl(:,ai);
#### 1.1.4.3. particle
xn_traj = zeros(nNonLin, N_P, N_T); % Collection of all trajectories
#### 1.1.4.4. how to use GP to predict
#### 1.1.4.5. sequential process
![alt text](figs/sequential.png)
# 2. Q&A
## 2.1. why does dynModel contain noise when predicting state?
At first I did not think it is good to use dynModel for propogate the state of particles, since it add some noise to the state, being so strange. However, now I accept that since the noise addition is how the pf gets its covariance or uncertainty.
## 2.2. SS = dy* P(:,:,i)*dy' + R;
dy is the derivative of measurement with respect to the map. The map containes 515 coefficients, including linear kernel and RBF kernel.
![alt text](figs/coeff.jpg)


# 3. Issue
# 4. Todo
## 4.1. compare smoother with filter

## 4.2. try to understand the dynModel: does it contain anything related to angular velocity?

## 4.3. 3.3
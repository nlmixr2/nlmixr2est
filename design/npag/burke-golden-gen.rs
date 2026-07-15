// Standalone golden-fixture generator for the Burke interior-point method.
// The `burke` function body is copied VERBATIM from pmcore 0.25.2
// (src/routines/estimation/ipm.rs) with only the Psi/Weights wrappers stubbed,
// so the numerics are identical to the Pmetrics reference.  Emits, for each test
// matrix, the input psi and the resulting weights + objective as CSV to stdout.
use faer::linalg::triangular_solve::solve_lower_triangular_in_place;
use faer::linalg::triangular_solve::solve_upper_triangular_in_place;
use faer::{Col, Mat, Row};
use anyhow::bail;
#[allow(unused_imports)]
use rayon::prelude::*;

pub struct Psi {
    matrix: Mat<f64>,
}
impl Psi {
    pub fn matrix(&self) -> &Mat<f64> {
        &self.matrix
    }
}
impl From<Mat<f64>> for Psi {
    fn from(matrix: Mat<f64>) -> Self {
        Psi { matrix }
    }
}

pub struct Weights {
    col: Col<f64>,
}
impl From<Col<f64>> for Weights {
    fn from(col: Col<f64>) -> Self {
        Weights { col }
    }
}
impl Weights {
    pub fn as_slice(&self) -> Vec<f64> {
        (0..self.col.nrows()).map(|i| *self.col.get(i)).collect()
    }
}

pub fn burke(psi: &Psi) -> anyhow::Result<(Weights, f64)> {
    let mut psi = psi.matrix().to_owned();

    psi.row_iter_mut().try_for_each(|row| {
        row.iter_mut().try_for_each(|x| {
            if !x.is_finite() {
                bail!("Input matrix must have finite entries")
            } else {
                *x = x.abs();
                Ok(())
            }
        })
    })?;

    let (n_sub, n_point) = psi.shape();

    let ecol: Col<f64> = Col::from_fn(n_point, |_| 1.0);
    let erow: Row<f64> = Row::from_fn(n_sub, |_| 1.0);

    let mut plam: Col<f64> = &psi * &ecol;
    let eps: f64 = 1e-8;
    let mut sig: f64 = 0.0;

    let mut lam = ecol.clone();

    let mut w: Col<f64> = Col::from_fn(plam.nrows(), |i| 1.0 / plam.get(i));

    let mut ptw: Col<f64> = psi.transpose() * &w;

    let ptw_max = ptw.iter().fold(f64::NEG_INFINITY, |acc, &x| x.max(acc));
    let shrink = 2.0 * ptw_max;
    lam *= shrink;
    plam *= shrink;
    w /= shrink;
    ptw /= shrink;

    let mut y: Col<f64> = &ecol - &ptw;
    let mut r: Col<f64> = Col::from_fn(n_sub, |i| erow.get(i) - w.get(i) * plam.get(i));
    let mut norm_r: f64 = r.iter().fold(0.0, |max, &val| max.max(val.abs()));

    let sum_log_plam: f64 = plam.iter().map(|x| x.ln()).sum();
    let sum_log_w: f64 = w.iter().map(|x| x.ln()).sum();
    let mut gap: f64 = (sum_log_w + sum_log_plam).abs() / (1.0 + sum_log_plam);

    let mut mu = lam.transpose() * &y / n_point as f64;

    let mut psi_inner: Mat<f64> = Mat::zeros(psi.nrows(), psi.ncols());
    let n_threads = faer::get_global_parallelism().degree();
    let rows = psi.nrows();
    let mut output: Vec<Mat<f64>> = (0..n_threads).map(|_| Mat::zeros(rows, rows)).collect();
    let mut h: Mat<f64> = Mat::zeros(rows, rows);

    while mu > eps || norm_r > eps || gap > eps {
        let smu = sig * mu;
        let inner = Col::from_fn(lam.nrows(), |i| lam.get(i) / y.get(i));
        let w_plam = Col::from_fn(plam.nrows(), |i| plam.get(i) / w.get(i));

        if psi.ncols() > n_threads * 128 {
            psi_inner
                .par_col_partition_mut(n_threads)
                .zip(psi.par_col_partition(n_threads))
                .zip(inner.par_partition(n_threads))
                .zip(output.par_iter_mut())
                .for_each(|(((mut psi_inner, psi), inner), output)| {
                    psi_inner
                        .as_mut()
                        .col_iter_mut()
                        .zip(psi.col_iter())
                        .zip(inner.iter())
                        .for_each(|((col, psi_col), inner_val)| {
                            col.iter_mut().zip(psi_col.iter()).for_each(|(x, psi_val)| {
                                *x = psi_val * inner_val;
                            });
                        });
                    faer::linalg::matmul::triangular::matmul(
                        output.as_mut(),
                        faer::linalg::matmul::triangular::BlockStructure::TriangularLower,
                        faer::Accum::Replace,
                        &psi_inner,
                        faer::linalg::matmul::triangular::BlockStructure::Rectangular,
                        psi.transpose(),
                        faer::linalg::matmul::triangular::BlockStructure::Rectangular,
                        1.0,
                        faer::Par::Seq,
                    );
                });
            let mut first_iter = true;
            for output in &output {
                if first_iter {
                    h.copy_from(output);
                    first_iter = false;
                } else {
                    h += output;
                }
            }
        } else {
            psi_inner
                .as_mut()
                .col_iter_mut()
                .zip(psi.col_iter())
                .zip(inner.iter())
                .for_each(|((col, psi_col), inner_val)| {
                    col.iter_mut().zip(psi_col.iter()).for_each(|(x, psi_val)| {
                        *x = psi_val * inner_val;
                    });
                });
            faer::linalg::matmul::triangular::matmul(
                h.as_mut(),
                faer::linalg::matmul::triangular::BlockStructure::TriangularLower,
                faer::Accum::Replace,
                &psi_inner,
                faer::linalg::matmul::triangular::BlockStructure::Rectangular,
                psi.transpose(),
                faer::linalg::matmul::triangular::BlockStructure::Rectangular,
                1.0,
                faer::Par::Seq,
            );
        }

        for i in 0..h.nrows() {
            h[(i, i)] += w_plam[i];
        }

        let uph = match h.llt(faer::Side::Lower) {
            Ok(llt) => llt,
            Err(_) => bail!("Error during Cholesky decomposition."),
        };
        let uph = uph.L().transpose().to_owned();

        let smuyinv: Col<f64> = Col::from_fn(ecol.nrows(), |i| smu * (ecol[i] / y[i]));
        let psi_dot_muyinv: Col<f64> = &psi * &smuyinv;
        let rhsdw: Row<f64> = Row::from_fn(erow.ncols(), |i| erow[i] / w[i] - psi_dot_muyinv[i]);
        let mut dw = Mat::from_fn(rhsdw.ncols(), 1, |i, _j| *rhsdw.get(i));

        solve_lower_triangular_in_place(uph.transpose().as_ref(), dw.as_mut(), faer::Par::rayon(0));
        solve_upper_triangular_in_place(uph.as_ref(), dw.as_mut(), faer::Par::rayon(0));

        let dw = dw.col(0);
        let dy = -(psi.transpose() * dw);
        let inner_times_dy = Col::from_fn(ecol.nrows(), |i| inner[i] * dy[i]);
        let dlam: Row<f64> = Row::from_fn(ecol.nrows(), |i| smuyinv[i] - lam[i] - inner_times_dy[i]);

        let ratio_dlam_lam = Row::from_fn(lam.nrows(), |i| dlam[i] / lam[i]);
        let min_ratio_dlam = ratio_dlam_lam.iter().cloned().fold(f64::INFINITY, f64::min);
        let mut alfpri: f64 = -1.0 / min_ratio_dlam.min(-0.5);
        alfpri = (0.99995 * alfpri).min(1.0);

        let ratio_dy_y = Row::from_fn(y.nrows(), |i| dy[i] / y[i]);
        let min_ratio_dy = ratio_dy_y.iter().cloned().fold(f64::INFINITY, f64::min);
        let ratio_dw_w = Row::from_fn(dw.nrows(), |i| dw[i] / w[i]);
        let min_ratio_dw = ratio_dw_w.iter().cloned().fold(f64::INFINITY, f64::min);
        let mut alfdual = -1.0 / min_ratio_dy.min(-0.5);
        alfdual = alfdual.min(-1.0 / min_ratio_dw.min(-0.5));
        alfdual = (0.99995 * alfdual).min(1.0);

        lam += alfpri * dlam.transpose();
        w += alfdual * dw;
        y += alfdual * &dy;

        mu = lam.transpose() * &y / n_point as f64;
        plam = &psi * &lam;
        r = Col::from_fn(n_sub, |i| erow.get(i) - w.get(i) * plam.get(i));
        ptw -= alfdual * dy;

        norm_r = r.norm_max();
        let sum_log_plam: f64 = plam.iter().map(|x| x.ln()).sum();
        let sum_log_w: f64 = w.iter().map(|x| x.ln()).sum();
        gap = (sum_log_w + sum_log_plam).abs() / (1.0 + sum_log_plam);

        if mu < eps && norm_r > eps {
            sig = 1.0;
        } else {
            let candidate1 = (1.0 - alfpri).powi(2);
            let candidate2 = (1.0 - alfdual).powi(2);
            let candidate3 = (norm_r - mu) / (norm_r + 100.0 * mu);
            sig = candidate1.max(candidate2).max(candidate3).min(0.3);
        }
    }
    lam /= n_sub as f64;
    let obj = (&psi * &lam).iter().map(|x| x.ln()).sum();
    let lam_sum: f64 = lam.iter().sum();
    lam = &lam / lam_sum;

    Ok((lam.into(), obj))
}

fn emit(name: &str, mat: &Mat<f64>) {
    let psi = Psi::from(mat.clone());
    let (w, obj) = burke(&psi).unwrap();
    let (nr, nc) = mat.shape();
    // input matrix
    for i in 0..nr {
        let row: Vec<String> = (0..nc).map(|j| format!("{:.17e}", mat[(i, j)])).collect();
        println!("{},psi,{},{}", name, i, row.join(","));
    }
    // weights
    let ws: Vec<String> = w.as_slice().iter().map(|x| format!("{:.17e}", x)).collect();
    println!("{},weights,-1,{}", name, ws.join(","));
    // objective
    println!("{},objective,-1,{:.17e}", name, obj);
}

fn main() {
    // Case 1: identity 5x5 -> uniform 1/5
    emit("identity5", &Mat::identity(5, 5));
    // Case 2: uniform 4x6 (all ones) -> uniform 1/6
    emit("uniform4x6", &Mat::from_fn(4, 6, |_, _| 1.0));
    // Case 3: non-uniform 4x5, deterministic entries
    emit("nonuniform4x5", &Mat::from_fn(4, 5, |i, j| {
        1.0 + (i as f64) * 0.5 + (j as f64) * (j as f64) * 0.3 + if j == 0 { 4.0 } else { 0.0 }
    }));
    // Case 4: realistic-ish 6x8 with a couple of dominant columns
    emit("mix6x8", &Mat::from_fn(6, 8, |i, j| {
        let base = 0.2 + 0.05 * (i as f64) + 0.01 * (j as f64);
        if j == 2 { base + 3.0 } else if j == 5 { base + 1.5 } else { base }
    }));
}

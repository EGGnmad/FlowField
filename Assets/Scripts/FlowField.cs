using UnityEngine;

namespace NeroAdventure.Physics.Field
{
    public class FlowField
    {
        private int length;
        public float diffuseRate;
        public float viscosityRate;

        private float[] density;
        private float[] density_prev;

        private float[] Vx;
        private float[] Vy;

        private float[] Vx_prev;
        private float[] Vy_prev;

        public FlowField(int length, float viscosityRate, float diffuseRate)
        {
            this.length = length;
            this.diffuseRate = diffuseRate;
            this.viscosityRate = viscosityRate;

            density = new float[length * length];
            density_prev = new float[length * length];

            Vx = new float[length * length];
            Vx_prev = new float[length * length];
            Vy = new float[length * length];
            Vy_prev = new float[length * length];
        }

        private int IX(int x, int y)
        {
            return Mathf.Clamp(x + y * length, 0, length * length - 1);
        }

        public void AddDensity(Vector2Int pos, float amount)
        {
            density[IX(pos.x, pos.y)] += amount;
        }

        public void AddVelocity(Vector2Int pos, Vector2 amount)
        {
            Vx[IX(pos.x, pos.y)] += amount.x;
            Vy[IX(pos.x, pos.y)] += amount.y;
        }

        private void Diffuse(int b, float[] x, float[] x0, float diff)
        {
            int N = length;
            float a = Time.deltaTime * diff * N * N;

            LinearSolve(b, x, x0, a, 1 + 4 * a);
        }

        private void Project(float[] velocX, float[] velocY, float[] p, float[] div)
        {
            int N = length;
            for (int j = 1; j < N - 1; j++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    div[IX(i, j)] = -0.5f * (
                        velocX[IX(i + 1, j)]
                        - velocX[IX(i - 1, j)]
                        + velocY[IX(i, j + 1)]
                        - velocY[IX(i, j - 1)]
                    ) / N;
                    p[IX(i, j)] = 0;
                }
            }

            SetBound(0, div);
            SetBound(0, p);
            LinearSolve(0, p, div, 1, 6);

            for (int j = 1; j < N - 1; j++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
                    velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
                }
            }

            SetBound(1, velocX);
            SetBound(2, velocY);
        }

        private void Advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY)
        {
            int N = length;
            float i0, i1, j0, j1;

            float dtx = Time.deltaTime * (N - 2);
            float dty = Time.deltaTime * (N - 2);

            float s0, s1, t0, t1;
            float tmp1, tmp2, x, y;

            float Nfloat = N;
            float ifloat, jfloat;
            int i, j;

            for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++)
            {
                for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++)
                {
                    tmp1 = dtx * velocX[IX(i, j)];
                    tmp2 = dty * velocY[IX(i, j)];
                    x = ifloat - tmp1;
                    y = jfloat - tmp2;

                    if (x < 0.5f) x = 0.5f;
                    if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
                    i0 = Mathf.Floor(x);
                    i1 = i0 + 1.0f;
                    if (y < 0.5f) y = 0.5f;
                    if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
                    j0 = Mathf.Floor(y);
                    j1 = j0 + 1.0f;

                    s1 = x - i0;
                    s0 = 1.0f - s1;
                    t1 = y - j0;
                    t0 = 1.0f - t1;

                    int i0i = (int)i0;
                    int i1i = (int)i1;
                    int j0i = (int)j0;
                    int j1i = (int)j1;

                    d[IX(i, j)] =
                        s0 * (t0 * d0[IX(i0i, j0i)] + (t1 * d0[IX(i0i, j1i)])) +
                        s1 * (t0 * d0[IX(i1i, j0i)] + (t1 * d0[IX(i1i, j1i)]));
                }
            }

            SetBound(b, d);
        }

        private void LinearSolve(int b, float[] x, float[] x0, float a, float c)
        {
            int N = length;

            float cRecip = 1.0f / c;
            for (int k = 0; k < 10; k++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    for (int i = 1; i < N - 1; i++)
                    {
                        x[IX(i, j)] =
                            (x0[IX(i, j)] +
                             a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip;
                    }
                }

                SetBound(b, x);
            }
        }

        private void SetBound(int b, float[] x)
        {
            int N = length;
            for (int i = 1; i < N - 1; i++)
            {
                x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
                x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
            }

            for (int j = 1; j < N - 1; j++)
            {
                x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
                x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
            }

            x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
            x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
            x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
            x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 1, N - 1)] + x[IX(N - 1, N - 2)]);
        }

        public void Step()
        {
            Diffuse(1, Vx_prev, Vx, viscosityRate);
            Diffuse(2, Vy_prev, Vy, viscosityRate);

            Project(Vx_prev, Vy_prev, Vx, Vy);

            Advect(1, Vx, Vx_prev, Vx_prev, Vy_prev);
            Advect(2, Vy, Vy_prev, Vx_prev, Vy_prev);

            Project(Vx, Vy, Vx_prev, Vy_prev);

            Diffuse(0, density_prev, density, diffuseRate);
            Advect(0, density, density_prev, Vx, Vy);
        }

        public int GetLength() => length;

        public float[,] GetDensityField()
        {
            float[,] densityField = new float[length, length];

            for (int j = 0; j < length; j++)
            {
                for (int i = 0; i < length; i++)
                {
                    densityField[j, i] = density[IX(i, j)];
                }
            }

            return densityField;
        }

        public void GetDensityField(ref float[,] field)
        {
            for (int j = 0; j < length; j++)
            {
                for (int i = 0; i < length; i++)
                {
                    field[j, i] = density[IX(i, j)];
                }
            }
        }
        
        public Vector2[,] GetVelocityField()
        {
            Vector2[,] velocityField = new Vector2[length, length];

            for (int j = 0; j < length; j++)
            {
                for (int i = 0; i < length; i++)
                {
                    velocityField[j, i] = new Vector2(Vx[IX(i, j)], Vy[IX(i, j)]);
                }
            }

            return velocityField;
        }

        // Better Performance
        public void GetVelocityField(ref Vector2[,] field)
        {
            for (int j = 0; j < length; j++)
            {
                for (int i = 0; i < length; i++)
                {
                    field[j, i] = new Vector2(Vx[IX(i, j)], Vy[IX(i, j)]);
                }
            }
        }
    }
}
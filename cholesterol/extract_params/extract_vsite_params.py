import random
import math
import numpy as np
import torch
from torch import nn


def read_data():
    data = []
    with open("cholesterol.dat") as infile:
        for line in infile:
            cols = line.split()
            x = float(cols[3])
            y = float(cols[4])
            z = float(cols[5])
            data.append([x, y, z])

    data = np.array(data)
    return np.reshape(data, newshape=(-1, 8, 3))


def prepare_data(data, to_build, i, j, k):
    y = data[:, to_build, :]
    ri = data[:, i, :]
    rij = data[:, j, :] - ri
    rik = data[:, k, :] - ri
    cross = np.cross(rij, rik)
    x = np.stack([ri, rij, rik, cross], axis=1)
    return torch.FloatTensor(x), torch.FloatTensor(y)


class VSite3(nn.Module):
    def __init__(self):
        super().__init__()

        self.a = nn.Parameter(torch.tensor(random.random() * 2 - 1))
        self.b = nn.Parameter(torch.tensor(random.random() * 2 - 1))
        self.c = nn.Parameter(torch.tensor(random.random() * 2 - 1))

    def forward(self, v):
        ri = v[:, 0, :]
        rij = v[:, 1, :]
        rik = v[:, 2, :]
        cross = v[:, 3, :]
        return ri + self.a * rij + self.b * rik + self.c * cross


# Read in the data
data = read_data()

# prepare data
to_build = 0
i = 3
j = 2
k = 6
x, y = prepare_data(data, to_build, i, j, k)

vsite = VSite3()
opt = torch.optim.Adam(vsite.parameters(), lr=3e-5)

for i in range(100_000):
    yhat = vsite(x)
    mse = torch.sum((yhat - y) ** 2, dim=1)
    loss = torch.mean(mse)
    opt.zero_grad()
    loss.backward()
    opt.step()

print(f"RMSD after fitting: {math.sqrt(loss.item()):.3f}.")
print(f"a = {vsite.a.item():.4f}")
print(f"b = {vsite.b.item():.4f}")
print(f"c = {vsite.c.item():.4f}")
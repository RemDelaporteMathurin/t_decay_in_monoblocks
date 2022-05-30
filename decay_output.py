import FESTIM as F
import fenics as f


class DecayXDMF(F.XDMFExport):
    def __init__(
        self, label="decay", filename=None, mode=1, checkpoint=True, folder=None
    ) -> None:
        super().__init__("retention", label, filename, mode, checkpoint, folder)

    def write(self, t):
        functionspace = self.function.function_space()
        density_as_function = f.project(self.make_decay(), functionspace)
        self.function = density_as_function
        super().write(t)

    def make_decay(self):
        tritium_half_life = 12.4  # years
        tritium_half_life *= 364.25  # days
        tritium_half_life *= 24  # hours
        tritium_half_life *= 3600  # seconds

        decay_constant = 0.69 / tritium_half_life  # ln(2)/half_life  in s-1
        decay_mobile = decay_constant * self.function

        return decay_mobile

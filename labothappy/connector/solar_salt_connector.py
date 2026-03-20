# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 09:11:54 2026

@author: gregoire.hendrix

Solar Salt (60% NaNO3 / 40% KNO3) connector.
Inherits from MassConnector and bypasses CoolProp entirely.
All thermophysical properties are computed from Zavoico (2001) correlations.

Reference:
    Zavoico, A.B. (2001). Solar Power Tower Design Basis Document.
    Sandia National Laboratories, SAND2001-2100.
    Valid range: 533 K (260°C) to 873 K (600°C).
"""

import math
import warnings
from labothappy.connector.mass_connector import MassConnector


class SolarSaltConnector(MassConnector):
    """
    MassConnector for Solar Salt (60% NaNO3 / 40% KNO3).

    CoolProp is bypassed entirely. All properties are derived from the
    Zavoico (2001) correlations using T as the primary state variable.

    Supported input pairs for calculate_properties:
        (T, P)  — temperature + pressure  [most common]
        (H, P)  — enthalpy  + pressure    [used by the recursive solver]

    Reference enthalpy: T_ref = 260°C (533.15 K), h_ref = 0 J/kg.
    """

    T_REF_K = 533.15   # reference temperature for enthalpy [K]

    def __init__(self):
        # Initialise parent without a fluid to prevent any CoolProp call
        super().__init__(fluid=None)
        # Label only — never passed to CoolProp
        self.fluid = 'SolarSalt'

    # ------------------------------------------------------------------
    # Zavoico (2001) static correlations
    # ------------------------------------------------------------------

    @staticmethod
    def _Cp(T_K):
        """
        Specific heat capacity [J/(kg·K)].
        Cp = 1443 + 0.172 * T_C
        """
        T_C = T_K - 273.15
        return 1443.0 + 0.172 * T_C

    @staticmethod
    def _rho(T_K):
        """
        Density [kg/m³].
        rho = 2090 - 0.636 * T_C
        """
        T_C = T_K - 273.15
        return 2090.0 - 0.636 * T_C

    @staticmethod
    def _lam(T_K):
        """
        Thermal conductivity [W/(m·K)].
        lambda = 0.443 + 1.9e-4 * T_C
        """
        T_C = T_K - 273.15
        return 0.443 + 1.9e-4 * T_C

    @staticmethod
    def _mu(T_K):
        """
        Dynamic viscosity [Pa·s].
        mu = 1e-3 * (22.714 - 0.120*T_C + 2.281e-4*T_C² - 1.474e-7*T_C³)
        """
        T_C = T_K - 273.15
        return 1e-3 * (22.714
                       - 0.120    * T_C
                       + 2.281e-4 * T_C**2
                       - 1.474e-7 * T_C**3)

    @classmethod
    def _h(cls, T_K):
        """
        Specific enthalpy [J/kg] relative to T_ref = 260°C (533.15 K).
        h = integral(Cp dT) from T_ref to T
          = 1443*(T_C - T_ref_C) + 0.172/2*(T_C² - T_ref_C²)
        """
        T_C     = T_K         - 273.15
        T_ref_C = cls.T_REF_K - 273.15
        return (1443.0 * (T_C - T_ref_C)
                + 0.172 / 2.0 * (T_C**2 - T_ref_C**2))

    @classmethod
    def _T_from_h(cls, h_Jkg):
        """
        Invert h(T) analytically to recover T [K] from h [J/kg].

        h = 1443*(T_C - T_ref_C) + 0.086*(T_C² - T_ref_C²)
        Rearranged as quadratic in T_C:
            0.086*T_C² + 1443*T_C - (h + 1443*T_ref_C + 0.086*T_ref_C²) = 0
        """
        T_ref_C = cls.T_REF_K - 273.15
        h_ref   = 1443.0 * T_ref_C + 0.086 * T_ref_C**2

        a =  0.086
        b =  1443.0
        c = -(h_Jkg + h_ref)

        discriminant = b**2 - 4.0 * a * c
        if discriminant < 0:
            raise ValueError(
                f"SolarSaltConnector._T_from_h: negative discriminant "
                f"for h = {h_Jkg:.2f} J/kg"
            )
        T_C = (-b + math.sqrt(discriminant)) / (2.0 * a)
        return T_C + 273.15   # [K]

    # ------------------------------------------------------------------
    # Override calculate_properties — no CoolProp involved
    # ------------------------------------------------------------------

    def calculate_properties(self):
        """
        Compute all thermophysical properties from the known state variables
        using Zavoico (2001) correlations.

        Supported input pairs: (T, P) or (H, P).
        """
        keys = [v[0] for v in self.variables_input]

        try:
            # --- 1. Resolve temperature ---
            if 'T' in keys:
                T_K = self.T
            elif 'H' in keys:
                T_K = self._T_from_h(self.h)
                self.T = T_K
            else:
                return   # insufficient information — wait for more inputs

            # --- 2. Compute all properties ---
            self.T  = T_K
            self.cp = self._Cp(T_K)
            self.D  = self._rho(T_K)
            self.h  = self._h(T_K)
            self.s  = None   # entropy not required by the cycle solver
            self.x  = None   # no phase change for molten salts

            # Pressure — keep user-set value or default to 1 bar
            if self.p is None:
                self.p = 1e5

            self.state_known = True

        except Exception as e:
            warnings.warn(f"SolarSaltConnector.calculate_properties failed: {e}")

    # ------------------------------------------------------------------
    # Override check_completely_known — bypass the fluid != None guard
    # ------------------------------------------------------------------

    def check_completely_known(self):
        """
        Same logic as MassConnector.check_completely_known but without the
        'fluid != None' guard that would short-circuit the calculation.
        """
        if len(self.variables_input) >= 2 or (
                len(self.variables_input) == 1
                and (self.SC is not None or self.SH is not None)):
            self.calculate_properties()

        if (self.m_dot is not None or self.V_dot is not None) and self.state_known:
            if self.m_dot is not None and self.D is not None:
                self.V_dot = self.m_dot / self.D * 3600
            self.completely_known = True
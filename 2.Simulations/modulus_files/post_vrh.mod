# Output final values

variable C11all equal ${C11}
variable C22all equal ${C22}
variable C33all equal ${C33}

variable C12all equal 0.5*(${C12}+${C21})
variable C13all equal 0.5*(${C13}+${C31})
variable C23all equal 0.5*(${C23}+${C32})

variable C44all equal ${C44}
variable C55all equal ${C55}
variable C66all equal ${C66}

variable C14all equal 0.5*(${C14}+${C41})
variable C15all equal 0.5*(${C15}+${C51})
variable C16all equal 0.5*(${C16}+${C61})

variable C24all equal 0.5*(${C24}+${C42})
variable C25all equal 0.5*(${C25}+${C52})
variable C26all equal 0.5*(${C26}+${C62})

variable C34all equal 0.5*(${C34}+${C43})
variable C35all equal 0.5*(${C35}+${C53})
variable C36all equal 0.5*(${C36}+${C63})

variable C45all equal 0.5*(${C45}+${C54})
variable C46all equal 0.5*(${C46}+${C64})
variable C56all equal 0.5*(${C56}+${C65})

# Cubic averages

variable C11cubic equal (${C11all}+${C22all}+${C33all})/3.0
variable C12cubic equal (${C12all}+${C13all}+${C23all})/3.0
variable C44cubic equal (${C44all}+${C55all}+${C66all})/3.0

# Isotropic-equivalent moduli (Voigt K & Voigt G)
variable bulkmodulus equal (${C11cubic}+2.0*${C12cubic})/3.0
variable shearmodulus1 equal ${C44cubic}
variable shearmodulus2 equal (${C11cubic}-${C12cubic})/2.0
variable shearmodulus equal ((${C11cubic}-${C12cubic})+3.0*${C44cubic})/5.0

# variable poissonratio equal 1.0/(1.0+${C11cubic}/${C12cubic})
# variable youngsmodulus equal 3.0*${bulkmodulus}*(1.0-2.0*${poissonratio})

# Poisson ratio and Young's modulus (isotropic relations)
variable poissonratio equal (3.0*${bulkmodulus}-2.0*${shearmodulus})/(2.0*(3.0*${bulkmodulus}+${shearmodulus}))
variable youngsmodulus equal (9.0*${bulkmodulus}*${shearmodulus})/(3.0*${bulkmodulus}+${shearmodulus})

print "========================================="
print "Components of the Elastic Constant Tensor"
print "========================================="

print "Elastic Constant C11all = ${C11all} ${cunits}"
print "Elastic Constant C22all = ${C22all} ${cunits}"
print "Elastic Constant C33all = ${C33all} ${cunits}"

print "Elastic Constant C12all = ${C12all} ${cunits}"
print "Elastic Constant C13all = ${C13all} ${cunits}"
print "Elastic Constant C23all = ${C23all} ${cunits}"

print "Elastic Constant C44all = ${C44all} ${cunits}"
print "Elastic Constant C55all = ${C55all} ${cunits}"
print "Elastic Constant C66all = ${C66all} ${cunits}"

print "Elastic Constant C14all = ${C14all} ${cunits}"
print "Elastic Constant C15all = ${C15all} ${cunits}"
print "Elastic Constant C16all = ${C16all} ${cunits}"

print "Elastic Constant C24all = ${C24all} ${cunits}"
print "Elastic Constant C25all = ${C25all} ${cunits}"
print "Elastic Constant C26all = ${C26all} ${cunits}"

print "Elastic Constant C34all = ${C34all} ${cunits}"
print "Elastic Constant C35all = ${C35all} ${cunits}"
print "Elastic Constant C36all = ${C36all} ${cunits}"

print "Elastic Constant C45all = ${C45all} ${cunits}"
print "Elastic Constant C46all = ${C46all} ${cunits}"
print "Elastic Constant C56all = ${C56all} ${cunits}"

print "========================================="
print "Isotropic-equivalent moduli (Voigt-based)"
print "========================================="

print "Bulk Modulus = ${bulkmodulus} ${cunits}"
print "Shear Modulus 1 = ${shearmodulus1} ${cunits}"
print "Shear Modulus 2 = ${shearmodulus2} ${cunits}"
print "Shear Modulus = ${shearmodulus} ${cunits}"
print "Poisson Ratio = ${poissonratio}"
print "Youngs Modulus = ${youngsmodulus} ${cunits}"
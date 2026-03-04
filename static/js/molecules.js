// Simple SMILES parser for basic visualization
function parseSMILES(smiles) {
    const atoms = []
    const bonds = []
    const branchStack = []
    const ringClosures = {}
    let currentIndex = 0
    let currentAtomIndex = -1
    let pendingBondType = 1

    function atomFromChar(char) {
        if (char === 'C' || char === 'O' || char === 'N' || char === 'S' || char === 'P' || char === 'H') {
            return char
        }
        return null
    }

    // Parse SMILES notation (hydrogens are implicit, not shown)
    for (let i = 0; i < smiles.length; i++) {
        const char = smiles[i]
        const atomSymbol = atomFromChar(char)

        if (atomSymbol) {
            // Keep hydrogens implicit in fallback view for readability
            if (atomSymbol !== 'H') {
                atoms.push({ element: atomSymbol, index: currentIndex })

                if (currentAtomIndex >= 0) {
                    bonds.push({ start: currentAtomIndex, end: currentIndex, type: pendingBondType })
                }

                currentAtomIndex = currentIndex
                currentIndex++
            }

        } else if (char === '(') {
            branchStack.push(currentAtomIndex)
        } else if (char === ')') {
            if (branchStack.length > 0) {
                currentAtomIndex = branchStack.pop()
            }
        } else if (char === '=') {
            pendingBondType = 2
        } else if (char === '#') {
            pendingBondType = 3
        } else if (char >= '0' && char <= '9') {
            if (currentAtomIndex >= 0) {
                if (ringClosures[char] === undefined) {
                    ringClosures[char] = {
                        atomIndex: currentAtomIndex,
                        bondType: pendingBondType,
                    }
                } else {
                    const start = ringClosures[char].atomIndex
                    const type = pendingBondType !== 1 ? pendingBondType : ringClosures[char].bondType
                    bonds.push({ start, end: currentAtomIndex, type })
                    delete ringClosures[char]
                }
            }
        }

        if (char !== '=') {
            if (char !== '#') {
                pendingBondType = 1
            }
        }
    }

    return { atoms, bonds }
}

// Create 3D molecule visualization
function createMoleculeViewer(containerId, smiles, atomsData, bondsData) {
    const container = document.getElementById(containerId)
    if (!container) return

    const width = container.clientWidth
    const height = 300

    // Setup scene
    const scene = new THREE.Scene()
    scene.background = new THREE.Color(0xffffff)

    // Setup camera
    const camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000)
    camera.position.z = 10

    // Setup renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true })
    renderer.setSize(width, height)
    container.appendChild(renderer.domElement)

    // Setup controls
    const controls = new THREE.OrbitControls(camera, renderer.domElement)
    controls.enableDamping = true
    controls.dampingFactor = 0.05

    // Add lights
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.5)
    scene.add(ambientLight)

    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8)
    directionalLight.position.set(5, 5, 5)
    scene.add(directionalLight)

    // Use backend-provided 3D coordinates or fallback to SMILES parsing
    let atoms, bonds, positions

    if (atomsData && atomsData.length > 0) {
        // Use backend-provided data with real 3D coordinates
        atoms = atomsData.map((atom, i) => ({ element: atom.element, index: i }))
        positions = atomsData.map(atom => ({ x: atom.x, y: atom.y, z: atom.z }))
        bonds = bondsData || []
    } else {
        // Fallback to SMILES parsing (legacy)
        const parsed = parseSMILES(smiles)
        atoms = parsed.atoms
        bonds = parsed.bonds

        // Generate simple 3D positions with generous spacing
        const angleStep = (2 * Math.PI) / Math.max(atoms.length, 1)
        const radius = Math.max(12, atoms.length * 2.5)

        positions = atoms.map((atom, i) => {
            const angle = i * angleStep
            return {
                x: radius * Math.cos(angle),
                y: radius * Math.sin(angle) * 0.8,
                z: (Math.random() - 0.5) * 2
            }
        })
    }

    // Create atoms
    atoms.forEach((atom, i) => {
        const geometry = new THREE.SphereGeometry(
            atomRadii[atom.element] || 0.5,
            32,
            32
        )
        const material = new THREE.MeshPhongMaterial({
            color: atomColors[atom.element] || 0x808080,
            shininess: 30
        })
        const sphere = new THREE.Mesh(geometry, material)
        sphere.position.set(positions[i].x, positions[i].y, positions[i].z)
        scene.add(sphere)

        // Add atom label (skip for hydrogen to reduce clutter)
        if (atom.element !== 'H') {
            // const canvas = document.createElement('canvas')
            // const context = canvas.getContext('2d')
            // canvas.width = 64
            // canvas.height = 64
            // context.fillStyle = '#000000'
            // context.font = 'Bold 48px Arial'
            // context.textAlign = 'center'
            // context.textBaseline = 'middle'
            // context.fillText(atom.element, 32, 32)

            // const texture = new THREE.CanvasTexture(canvas)
            // const spriteMaterial = new THREE.SpriteMaterial({ map: texture })
            // const sprite = new THREE.Sprite(spriteMaterial)
            // sprite.position.set(positions[i].x, positions[i].y + 1, positions[i].z)
            // sprite.scale.set(1, 1, 1)
            // scene.add(sprite)
        }
    })

    function createBondCylinder(startVec, endVec, radius, color) {
        const direction = new THREE.Vector3().subVectors(endVec, startVec)
        const length = direction.length()
        if (length < 1e-6) return

        const geometry = new THREE.CylinderGeometry(radius, radius, length, 14)
        const material = new THREE.MeshPhongMaterial({ color, shininess: 20 })
        const cylinder = new THREE.Mesh(geometry, material)

        const midpoint = new THREE.Vector3().addVectors(startVec, endVec).multiplyScalar(0.5)
        cylinder.position.copy(midpoint)

        const yAxis = new THREE.Vector3(0, 1, 0)
        cylinder.quaternion.setFromUnitVectors(yAxis, direction.clone().normalize())
        scene.add(cylinder)
    }

    function drawBondSet(startPos, endPos, bondType) {
        const direction = new THREE.Vector3().subVectors(endPos, startPos)
        const length = direction.length()
        if (length < 1e-6) return

        const unitDirection = direction.clone().normalize()
        const up = new THREE.Vector3(0, 1, 0)
        const alternate = new THREE.Vector3(1, 0, 0)
        const refAxis = Math.abs(unitDirection.dot(up)) > 0.95 ? alternate : up
        const perp = new THREE.Vector3().crossVectors(unitDirection, refAxis).normalize()

        if (bondType >= 2.5) {
            const offsets = [-0.18, 0, 0.18]
            offsets.forEach((offset) => {
                const shift = perp.clone().multiplyScalar(offset)
                createBondCylinder(
                    startPos.clone().add(shift),
                    endPos.clone().add(shift),
                    0.05,
                    0xFF6600
                )
            })
        } else if (bondType >= 1.5) {
            const offsets = [-0.12, 0.12]
            offsets.forEach((offset) => {
                const shift = perp.clone().multiplyScalar(offset)
                createBondCylinder(
                    startPos.clone().add(shift),
                    endPos.clone().add(shift),
                    0.045,
                    0x000000
                )
            })
        } else {
            createBondCylinder(startPos, endPos, 0.05, 0x000000)
        }
    }

    // Create bonds
    bonds.forEach((bond) => {
        // Handle both array format [start, end] and object format {start, end, type}
        const start = Array.isArray(bond) ? bond[0] : bond.start
        const end = Array.isArray(bond) ? bond[1] : bond.end
        const bondType = Number(Array.isArray(bond) ? 1 : (bond.type || 1))

        if (start < positions.length && end < positions.length) {
            const startPos = new THREE.Vector3(positions[start].x, positions[start].y, positions[start].z)
            const endPos = new THREE.Vector3(positions[end].x, positions[end].y, positions[end].z)
            drawBondSet(startPos, endPos, bondType)
        }
    })

    // Animation loop
    function animate() {
        requestAnimationFrame(animate)
        controls.update()
        renderer.render(scene, camera)
    }
    animate()

    // Handle resize
    window.addEventListener('resize', () => {
        const newWidth = container.clientWidth
        camera.aspect = newWidth / height
        camera.updateProjectionMatrix()
        renderer.setSize(newWidth, height)
    })
}

// Initialize all molecule viewers
window.addEventListener('load', function () {
    // Ensure THREE is loaded
    if (typeof THREE !== 'undefined') {
        initializeMolecules()
    } else {
        console.error('THREE.js failed to load')
    }
})

// Also initialize when loaded via HTMX
document.body.addEventListener('htmx:afterSettle', function (event) {
    if (typeof THREE !== 'undefined') {
        initializeMolecules()
    }
})

function initializeMolecules() {
    document.querySelectorAll('.molecule-viewer').forEach((container) => {
        container.innerHTML = ''
        const smiles = container.dataset.smiles
        const atomsData = container.dataset.atoms ? JSON.parse(container.dataset.atoms) : []
        const bondsData = container.dataset.bonds ? JSON.parse(container.dataset.bonds) : []
        createMoleculeViewer(container.id, smiles, atomsData, bondsData)
    })
}

// Expose functions globally for HTMX compatibility
window.initializeMolecules = initializeMolecules

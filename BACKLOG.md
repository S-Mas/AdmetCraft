# Backlog `v0.3.0` (Estado de mejoras post-`v0.2.1`)

Seguimiento de los 14 hallazgos no bloqueantes registrados tras `v0.2.0`/`v0.2.1`.

## 1. Conformidad JSON-RPC / MCP
- [x] **1.** Devolver error JSON-RPC `-32602` (*Invalid params*) ante parámetros inválidos en lugar de respuestas de herramienta genéricas (`apps/admetcraft-cli/src/mcp.rs`).
- [x] **2.** Validar estrictamente `"jsonrpc": "2.0"` y rechazar peticiones sin `jsonrpc` o con `"1.0"` con `-32600` (`apps/admetcraft-cli/src/mcp.rs`).
- [x] **3.** Distinguir notificaciones sin `id` (sin respuesta) frente a `"id": null` explícito (eco de `"id": null` según JSON-RPC 2.0).
- [x] **4.** Exigir `format` explícitamente en `export_report` (sin *fallback* silencioso a `csv`) y declarar `outputSchema` en todas las herramientas MCP (ahora 5 herramientas incluyendo `compare_molecules`).

## 2. CLI y ergonomía de sandbox
- [x] **5.** Validar `max_lipinski_violations` (y `min_qed` / `min_cns_mpo` / booleanos) rechazando valores negativos como `-1` con `-32602`.
- [x] **6.** Rechazar `--write-root` inexistente por defecto y exigir `--create-write-root` explícito (`PathSandbox::new_with_options`).
- [x] **7.** Emitir aviso explícito en `stderr` en modo `batch` sin `--read-root`/`--write-root`.

## 3. Informes y GUI
- [x] **8.** Truncado de etiquetas largas (`36` caracteres) en el SVG BOILED-Egg (`crates/admet-report/src/svg.rs`).
- [x] **9.** Resincronización del sketcher 2D tras `FilterBatch`, `ClearSession`, `LoadExampleBatch`, `ImportSdf`, `ImportCsv`, `SelectMolecule`, `Undo` y `Redo` (`apps/admetcraft/src/main.rs`).
- [x] **10.** Llamada a `ctx.request_repaint()` desde el hilo de exportación asíncrona (`apps/admetcraft/src/main.rs`) + selector nativo de archivos (`zenity`/`kdialog`).
- [ ] **11.** Soportar historial *Undo/Redo* unificado para edición carácter a carácter del cuadro multilínea SMILES cuando se alterna con comandos de sesión (diferido a `v0.4.0`).
- [ ] **12.** Unificar el renderizador vectorial 2D entre GUI (`egui`), SVG y PDF (diferido a `v0.4.0`).
- [x] **13.** Comprimir streams PDF con zlib `/Filter /FlateDecode` (~4.94× reducción de tamaño verificada con `zlib.decompress`).
- [x] **14.** Filtrar *noncharacters* XML (`U+FDD0..=U+FDEF`, `U+FFFE`, `U+FFFF`) además de los caracteres de control ASCII en `escape_xml`.

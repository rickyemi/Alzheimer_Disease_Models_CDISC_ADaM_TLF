# Screenshot the running Shiny app (start it first: shiny::runApp("app", port = 8765)). Requires: pip install playwright
import asyncio, sys
from playwright.async_api import async_playwright
import pathlib
OUT = str(pathlib.Path(__file__).resolve().parents[1] / "figures") + "/"
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        pg = await b.new_page(viewport={"width": 1500, "height": 1000}, device_scale_factor=1.5)
        errs = []
        pg.on("console", lambda m: errs.append(m.text) if m.type == "error" else None)
        await pg.goto("http://127.0.0.1:8765", wait_until="networkidle")
        await pg.wait_for_timeout(6000)
        await pg.screenshot(path=OUT + "22_dashboard_overview.png", full_page=True)
        async def tab(name, fname, wait=6000, setup=None):
            await pg.locator("a.nav-link", has_text=name).first.click()
            await pg.wait_for_timeout(1500)
            if setup: await setup()
            await pg.wait_for_timeout(wait)
            await pg.screenshot(path=OUT + fname, full_page=True)
        await tab("Central monitoring", "23_dashboard_central_monitoring.png")
        await tab("Safety", "_tmp1.png")
        await tab("Efficacy & biomarkers", "_tmp2.png")
        async def pick():
            await pg.get_by_label("Potential Hy's law").check()
            await pg.wait_for_timeout(2500)
        await tab("Patient profile", "24_dashboard_patient_profile.png", wait=9000, setup=pick)
        await tab("CDISC explorer", "_tmp3.png")
        txt = await pg.locator(".shiny-output-error").all_inner_texts()
        print("shiny errors:", txt[:5]); print("console errors:", errs[:5])
        await b.close()
asyncio.run(main())
